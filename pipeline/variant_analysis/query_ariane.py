"""
Classify every variant with the ARIANE service, storing results on the filesystem.

ARIANE (https://ariane-app.duckdns.org/) is an ENIGMA BRCA1/2 ACMG
classification service. This script reads gene + HGVS cDNA from the variant
table, POSTs them to the v1 batch endpoint, and writes one JSON file per variant.
Nothing is written to the database.

Output layout (under --out-dir):
    results/<VRS_Digest with ':' replaced by '_'>.json
    errors.tsv                variants the service could not classify: digest,
                              gene, c_notation, category, status, code,
                              retryable, message
    unqueryable.tsv           variants never sent, with the reason
    ariane_summary.tsv        one row per classified variant
    run.log                   (when launched detached)

Variants that already have a result file are skipped unless --overwrite, and so
are variants whose earlier failure the service marked non-retryable unless
--retry-failed. The run is idempotent, and a multi-day pass can be stopped and
resumed without re-spending quota on answers already known.

The API key is read from $ARIANE_API_KEY, else from --api-key-file (default
~/.config/ariane/api_key). It is deliberately not a command-line option -- argv
is visible to every user in the process list -- and it is never logged or
written to results.

Properties of the v1 service (1.9.11) that shape this script:

  * Quota is 5,000 classifications per key per UTC day, and failed
    classifications count against it. Every response reports what is left; when
    it runs out, all workers sleep until the reset rather than burning retries.
  * Errors come back per item inside a 200, each marked retryable or not, so a
    bad variant no longer takes its whole batch down. 422 bisection is kept as a
    safety net should a whole payload ever be rejected again.
  * A 504 does not mean the work was lost: the backend completes and caches it,
    so re-POSTing the same batch returns quickly.
  * A 401/403 means the key itself is missing or rejected, which retrying cannot
    fix, so the run stops at once rather than recording the same rejection for
    every remaining variant.
  * The deployed build revision is published only by /api/resources, not per
    result, so it is fetched at start, hourly, and whenever the per-result
    application version changes, and stamped on each result.
"""

import collections
import concurrent.futures
import csv
import datetime
import json
import logging
import os
import random
import re
import sys
import threading
import time

import click
import psycopg2
import requests

log = logging.getLogger(__name__)

DEFAULT_BASE_URL = 'https://ariane-app.duckdns.org/api/v1'
# The build revision is published here; v1 does not report it per result.
RESOURCES_URL = 'https://ariane-app.duckdns.org/api/resources'
DEFAULT_API_KEY_FILE = '~/.config/ariane/api_key'
API_KEY_HEADER = 'X-ARIANE-API-Key'

# The ARIANE developers recommend batches of 5 -- also published as
# recommended_uncached_batch_items in /api/v1/capabilities. 10 is the maximum.
DEFAULT_BATCH_SIZE = 5
MAX_BATCH_SIZE = 10
# concurrent_classifications_per_key in /api/v1/capabilities.
DEFAULT_WORKERS = 2

REQUEST_TIMEOUT = 90        # > the 60s gateway timeout, so we see the 504 itself
MAX_ATTEMPTS = 4
MAX_CONSECUTIVE_429 = 100   # a quota wait is not a failure, but never loop forever
INTER_REQUEST_DELAY = 0.5   # politeness: this is a small shared instance
BUILD_REFRESH_SECONDS = 3600
QUOTA_RESET_MARGIN = 30     # seconds past the advertised reset before resuming

_QUERY = """
SELECT "VRS_Digest", "Gene_Symbol", "HGVS_cDNA"
FROM variant
WHERE "Gene_Symbol" = ANY(%(genes)s)
"""

# Service-side evidence gaps. The service cannot classify these variants, and
# retrying will not change that.
_SERVICE_LIMITATION_CODES = {
    'structural_population_evidence_unavailable',
    'spliceai_coordinates_unavailable',
}


def utc_now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds')


def categorize_error(status, message, code=None, retryable=None):
    """Bucket a failure so the expected ones do not hide real breakage.

    service_limitation  the service cannot classify the variant: its reference
                        bundle cannot verify the reference allele, or evidence
                        it requires is unavailable
    transient           worth trying again on the next run
    data                the variant itself is unusable (malformed notation)
    error               anything else
    """
    text = str(message or '')
    if 'Reference allele could not be verified' in text or code in _SERVICE_LIMITATION_CODES:
        return 'service_limitation'
    # A SpliceAI lookup timing out is the Broad API being flaky under load, not a
    # property of the variant, whatever the item's retryable flag says.
    if retryable or status == 'retries_exhausted' or (
            'SpliceAI' in text and '(api_error)' in text):
        return 'transient'
    if code == 'invalid_variant' or status in ('skipped', 'missing'):
        return 'data'
    return 'error'


def normalize_notation(hgvs_cdna):
    """Reduce a stored HGVS cDNA string to the bare c. notation ARIANE wants.

    Some rows carry a transcript prefix ('NM_007294.3:c.*6207C>T'). Those are
    well-formed, just qualified, so strip the prefix and let the service judge
    them rather than discarding them as unqueryable. Returns None when there is
    no c. notation at all ('-' is the table's placeholder for missing).
    """
    if not hgvs_cdna:
        return None
    text = hgvs_cdna.strip()
    if ':c.' in text:
        text = text.split(':', 1)[1]
    return text if text.startswith('c.') else None


# ARIANE's installed reference bundle covers the coding sequence plus this many
# bases of intronic flank, and no UTR sequence at all. Outside it there is no
# reference base to verify, so a substitution is rejected outright while an
# indel still classifies (it needs no reference base). Measured in a live run:
# intronic substitutions at offsets 1-50 classified 3,468 of 3,468, at >=51
# failed 28,375 of 28,375; UTR substitutions failed 2,978 of 2,978.
DEFAULT_MAX_INTRON_OFFSET = 50

_INTRON_OFFSET_RE = re.compile(r'^c\.[0-9]+[+-]([0-9]+)')
_SUBSTITUTION_RE = re.compile(r'[ACGT]>[ACGT]$')


def intron_offset(c_notation):
    """Distance into the intron for a c.NNN+M / c.NNN-M notation, else None."""
    m = _INTRON_OFFSET_RE.match(c_notation or '')
    return int(m.group(1)) if m else None


def unverifiable_reason(c_notation, max_intron_offset=DEFAULT_MAX_INTRON_OFFSET,
                        skip_utr_substitutions=True):
    """Why the service cannot classify this variant, or None if it can.

    Sending a variant the reference bundle demonstrably cannot verify spends a
    unit of the daily quota on a guaranteed rejection, so those are filtered out
    before the run. Only substitutions are affected: indels carry no reference
    allele to check.
    """
    if not _SUBSTITUTION_RE.search(c_notation or ''):
        return None
    if skip_utr_substitutions and c_notation.startswith(('c.*', 'c.-')):
        return 'UTR substitution (outside ARIANE reference bundle)'
    if max_intron_offset is None:
        return None
    offset = intron_offset(c_notation)
    if offset is not None and offset > max_intron_offset:
        return (f'deep intronic substitution, offset {offset} > '
                f'{max_intron_offset} (outside ARIANE reference bundle)')
    return None


def digest_to_filename(vrs_digest):
    """Filesystem-safe filename for a VRS digest."""
    return vrs_digest.replace(':', '_') + '.json'


def load_api_key(key_file=DEFAULT_API_KEY_FILE):
    """The ARIANE API key: $ARIANE_API_KEY if set, else the key file.

    Deliberately not a command-line option -- argv is visible to every user in
    the process list. Errors name the file, never the key.
    """
    key = os.environ.get('ARIANE_API_KEY', '').strip()
    if key:
        return key
    path = os.path.expanduser(key_file)
    try:
        with open(path) as f:
            key = f.read().strip()
    except FileNotFoundError:
        raise click.ClickException(
            f'No ARIANE API key: set ARIANE_API_KEY or create {path} (mode 600).')
    if not key:
        raise click.ClickException(f'The ARIANE API key file {path} is empty.')
    mode = os.stat(path).st_mode & 0o777
    if mode & 0o077:
        log.warning('%s is readable by group or others (mode %o); chmod 600 it.',
                    path, mode)
    return key


class AuthError(Exception):
    """The API key is missing or rejected; retrying cannot help."""


class BuildTracker:
    """The deployed ARIANE build, which v1 reports only via /api/resources.

    Refreshed at first use, hourly, and whenever the per-result application
    version changes -- ARIANE has been redeployed mid-run before, and a single
    start-of-run check would mislabel everything after a redeploy. The build is
    therefore accurate to within the refresh interval; the exact per-result
    identity is the classifier_fingerprint ARIANE returns with every result.
    """

    def __init__(self, url=RESOURCES_URL, refresh_seconds=BUILD_REFRESH_SECONDS,
                 api_key=None):
        self.url = url
        self.refresh_seconds = refresh_seconds
        self.api_key = api_key
        self.build = None
        self.checked_at = None
        self._checked_mono = None
        self._seen_version = None
        self._lock = threading.Lock()

    def _fetch(self):
        # /api/resources answered without a key until 2026-09-11 and now
        # returns 401 without one, so the key is sent here too.
        headers = {API_KEY_HEADER: self.api_key} if self.api_key else {}
        try:
            r = requests.get(self.url, headers=headers, timeout=30)
            r.raise_for_status()
            return r.json().get('build_revision')
        except (requests.RequestException, ValueError) as e:
            log.warning('Could not read the ARIANE build revision from %s: %s',
                        self.url, e)
            return None

    def current(self, seen_version=None):
        """(build, checked_at) for results classified now."""
        with self._lock:
            version_changed = (seen_version is not None
                               and self._seen_version is not None
                               and seen_version != self._seen_version)
            if seen_version is not None:
                self._seen_version = seen_version
            if (self._checked_mono is None or version_changed
                    or time.monotonic() - self._checked_mono >= self.refresh_seconds):
                build = self._fetch()
                if build:
                    self.build = build
                self.checked_at = utc_now()
                self._checked_mono = time.monotonic()
            return self.build, self.checked_at


class Ariane:
    """Client for the ARIANE v1 batch classification endpoint."""

    def __init__(self, base_url, api_key=None, debug=False):
        self.base_url = base_url.rstrip('/')
        self.api_key = api_key
        self.debug = debug
        self._local = threading.local()
        # Shared by all workers: the daily quota the service last reported.
        self._quota_lock = threading.Lock()
        self.quota_remaining = None
        self.quota_reset = None          # epoch seconds

    @property
    def session(self):
        # One session per worker thread; requests.Session is not thread-safe.
        if not hasattr(self._local, 'session'):
            session = requests.Session()
            if self.api_key:
                session.headers[API_KEY_HEADER] = self.api_key
            self._local.session = session
        return self._local.session

    def _post(self, variants):
        payload = {'variants': [{'gene': g, 'c_notation': c} for _, g, c in variants]}
        return self.session.post(f'{self.base_url}/classify/batch',
                                 json=payload, timeout=REQUEST_TIMEOUT)

    def _note_quota(self, resp):
        """Record the daily quota the service reports on every response."""
        remaining = str(resp.headers.get('x-ratelimit-remaining', ''))
        reset = str(resp.headers.get('x-ratelimit-reset', ''))
        with self._quota_lock:
            if remaining.isdigit():
                self.quota_remaining = int(remaining)
            if reset.isdigit():
                self.quota_reset = int(reset)

    def wait_for_quota(self, needed):
        """Sleep until the reset if today's quota cannot cover `needed`."""
        with self._quota_lock:
            remaining, reset = self.quota_remaining, self.quota_reset
        if remaining is None or reset is None or remaining >= needed:
            return
        delay = reset - time.time() + QUOTA_RESET_MARGIN
        if delay > 0:
            log.warning('ARIANE daily quota nearly spent (%d left, %d needed); '
                        'sleeping %.1fh until the reset at 00:00 UTC.',
                        remaining, needed, delay / 3600)
            time.sleep(delay)
        with self._quota_lock:
            self.quota_remaining = None     # a new window; the next response says
                                            # how much of it is left

    def _wait_after_429(self, resp):
        """Wait out a 429: until the daily reset if the quota is spent, else as asked."""
        remaining = str(resp.headers.get('x-ratelimit-remaining', ''))
        reset = str(resp.headers.get('x-ratelimit-reset', ''))
        if remaining == '0' and reset.isdigit():
            delay = max(int(reset) - time.time(), 0) + QUOTA_RESET_MARGIN
            log.warning('ARIANE daily quota spent; sleeping %.1fh until the reset.',
                        delay / 3600)
        else:
            retry_after = str(resp.headers.get('Retry-After', ''))
            delay = float(retry_after) if retry_after.isdigit() else 60.0
            log.warning('ARIANE rate limit hit; sleeping %.0fs.', delay)
        time.sleep(delay)
        with self._quota_lock:
            self.quota_remaining = None

    def classify(self, variants, on_error):
        """Classify a list of (digest, gene, c_notation).

        Returns {digest: (classification, item_metadata)}. Variants the service
        rejects are passed to on_error(digest, gene, c_notation, status,
        message, code=..., retryable=...) instead. Raises AuthError if the key
        is rejected.
        """
        attempt = throttled = 0
        while attempt < MAX_ATTEMPTS and throttled < MAX_CONSECUTIVE_429:
            self.wait_for_quota(len(variants))
            try:
                resp = self._post(variants)
            except requests.RequestException as e:
                attempt += 1
                self._backoff(attempt, f'{type(e).__name__}: {e}')
                continue
            self._note_quota(resp)

            if resp.status_code == 429:
                # Out of quota or over the per-minute rate. Waiting it out is not
                # a failed attempt.
                throttled += 1
                self._wait_after_429(resp)
                continue
            throttled = 0

            if resp.status_code == 200:
                return self._collect(variants, resp.json(), on_error)

            if resp.status_code in (401, 403):
                raise AuthError(self._error_message(resp))

            if resp.status_code == 422:
                # The whole payload was rejected because of at least one bad
                # variant. Split to find it; the healthy ones are cached by now,
                # so the extra round trips are cheap.
                return self._bisect(variants, resp, on_error)

            attempt += 1
            if resp.status_code == 504:
                # The backend finished and cached the work even though nginx
                # gave up; retrying returns it quickly.
                log.info('504 on %d variant(s); retrying to collect cached results',
                         len(variants))
                self._backoff(attempt, '504 gateway timeout', base=2.0)
                continue

            if resp.status_code in (500, 502, 503):
                retry_after = resp.headers.get('Retry-After')
                delay = float(retry_after) if (retry_after or '').isdigit() else None
                self._backoff(attempt, f'HTTP {resp.status_code}', fixed=delay)
                continue

            # Anything else is not worth retrying.
            for digest, gene, c in variants:
                on_error(digest, gene, c, resp.status_code, self._error_message(resp))
            return {}

        for digest, gene, c in variants:
            on_error(digest, gene, c, 'retries_exhausted',
                     f'no success after {MAX_ATTEMPTS} attempts')
        return {}

    def _collect(self, variants, body, on_error):
        """Pull per-item results out of a 200 response."""
        results = {}
        by_index = {r.get('index'): r for r in body.get('results', [])}
        for i, (digest, gene, c) in enumerate(variants):
            item = by_index.get(i)
            if item is None:
                on_error(digest, gene, c, 'missing', 'no result at this index')
                continue
            # v1 returns the payload as "classification"; the pre-v1 endpoint
            # called it "result".
            payload = item.get('classification', item.get('result'))
            if item.get('status') == 'ok' and payload is not None:
                results[digest] = (payload, item.get('metadata') or {})
                continue
            err = item.get('error')
            if isinstance(err, dict):
                on_error(digest, gene, c, item.get('status', 'error'),
                         err.get('message', ''), code=err.get('code'),
                         retryable=err.get('retryable'))
            else:
                on_error(digest, gene, c, item.get('status', 'error'),
                         err or 'no result returned')
        return results

    def _bisect(self, variants, resp, on_error):
        """A 422 rejects the entire payload. Split until the offender is alone."""
        if len(variants) == 1:
            digest, gene, c = variants[0]
            on_error(digest, gene, c, 422, self._error_message(resp))
            return {}
        mid = len(variants) // 2
        out = {}
        for half in (variants[:mid], variants[mid:]):
            out.update(self.classify(half, on_error))
        return out

    @staticmethod
    def _error_message(resp):
        """The service's own explanation of a failed request."""
        try:
            body = resp.json()
        except ValueError:
            return resp.text[:300]
        if not isinstance(body, dict):
            return str(body)[:300]
        err = body.get('error')
        if isinstance(err, dict):
            return f"{err.get('code', '')}: {err.get('message', '')}"[:300]
        detail = body.get('detail')
        if isinstance(detail, list) and detail:
            return str(detail[0].get('msg', ''))[:300]
        return str(detail or body)[:300]

    def _backoff(self, attempt, reason, base=1.5, fixed=None):
        delay = fixed if fixed is not None else base * (2 ** (attempt - 1))
        delay += random.uniform(0, 0.5)     # jitter, so workers do not sync up
        log.warning('attempt %d failed (%s); sleeping %.1fs', attempt, reason, delay)
        time.sleep(delay)


def load_done(results_dir):
    """Digests already fetched -- the run's checkpoint."""
    return {f[:-len('.json')] for f in os.listdir(results_dir) if f.endswith('.json')}


def load_permanent_failures(errors_path):
    """Digests whose recorded failure the service marked non-retryable.

    Under a daily quota, resending them on every restart spends real
    classifications on answers already known. Transient failures are retried
    whatever their flag says, and rows from before errors.tsv carried a
    retryable column are ignored.
    """
    failed = set()
    if not os.path.exists(errors_path):
        return failed
    with open(errors_path) as f:
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 8 and parts[6] == 'false' and parts[3] != 'transient':
                failed.add(parts[0])
    return failed


def result_payload(digest, gene, c_notation, source, classification, item_meta,
                   build, build_checked_at):
    """The JSON written for one classified variant.

    _meta.ariane holds the per-result metadata ARIANE returned (application
    version, classifier fingerprint, policy, request id) plus the build revision
    current when the variant was classified.
    """
    ariane = dict(item_meta or {})
    ariane['build_version'] = build
    ariane['build_version_checked_at'] = build_checked_at
    return {
        '_meta': {'VRS_Digest': digest, 'gene': gene, 'c_notation': c_notation,
                  'source': source, 'retrieved': utc_now(), 'ariane': ariane},
        'result': classification,
    }


def write_summary(results_dir, summary_path):
    """One row per classified variant, from the stored JSON."""
    rows = []
    for name in sorted(os.listdir(results_dir)):
        if not name.endswith('.json'):
            continue
        with open(os.path.join(results_dir, name)) as f:
            d = json.load(f)
        meta, res = d.get('_meta', {}), d.get('result', {})
        ariane = meta.get('ariane') or {}
        # Each criterion is {name, applies, strength, points, reason, ...};
        # only the ones that actually applied are worth carrying here.
        applied = [c for c in (res.get('criteria') or []) if c.get('applies')]
        rows.append([
            meta.get('VRS_Digest'), meta.get('gene'), meta.get('c_notation'),
            res.get('predicted_class'), res.get('predicted_label'),
            res.get('total_points'), res.get('evidence_direction'),
            ';'.join(f'{c.get("name")}({c.get("strength")})' for c in applied),
            ariane.get('application_version'), ariane.get('build_version'),
        ])
    with open(summary_path, 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['VRS_Digest', 'gene', 'c_notation', 'predicted_class',
                    'predicted_label', 'total_points', 'evidence_direction',
                    'criteria_met', 'ariane_version', 'build_version'])
        w.writerows(rows)
    print(f'Summary: {len(rows)} variant(s) -> {summary_path}')


@click.command()
@click.option('--db-url', default='postgresql://postgres:postgres@localhost/storage.pg',
              envvar='PIPELINE_DB_URL', show_default=True)
@click.option('--schema', default='pipeline', show_default=True)
@click.option('--out-dir', default='/data/new_schema/data_out/ariane',
              show_default=True, type=click.Path(),
              help='Directory for results/, errors.tsv and the summary')
@click.option('--base-url', default=DEFAULT_BASE_URL, show_default=True)
@click.option('--api-key-file', default=DEFAULT_API_KEY_FILE, show_default=True,
              help='File holding the ARIANE API key (mode 600); $ARIANE_API_KEY '
                   'takes precedence. The key itself is never a command-line '
                   'argument.')
@click.option('--genes', default='BRCA1,BRCA2', show_default=True)
@click.option('--workers', default=DEFAULT_WORKERS, show_default=True,
              help='Concurrent requests. ARIANE allows 2 concurrent '
                   'classifications per key, and the daily quota, not '
                   'throughput, is what bounds the run.')
@click.option('--max-intron-offset', default=DEFAULT_MAX_INTRON_OFFSET,
              show_default=True, type=int,
              help='Skip intronic substitutions deeper than this many bases into '
                   'the intron: ARIANE\'s reference bundle cannot verify them, so '
                   'they are a guaranteed rejection. Pass -1 to send them anyway.')
@click.option('--send-utr-substitutions', is_flag=True, default=False,
              help='Send UTR substitutions anyway. By default they are skipped: '
                   'the reference bundle holds no UTR sequence, so each is a '
                   'guaranteed rejection that still costs quota.')
@click.option('--batch-size', default=DEFAULT_BATCH_SIZE, show_default=True,
              type=click.IntRange(1, MAX_BATCH_SIZE),
              help='Variants per request. The ARIANE developers recommend 5; the '
                   f'API accepts at most {MAX_BATCH_SIZE}.')
@click.option('--retry-failed', is_flag=True, default=False,
              help='Also resend variants whose earlier failure the service marked '
                   'non-retryable. Costs quota.')
@click.option('--vrs-digest', default=None, metavar='DIGEST',
              help='Classify only this variant')
@click.option('--limit', default=None, type=int, help='Classify at most this many')
@click.option('--overwrite', is_flag=True, default=False,
              help='Re-fetch variants that already have a result file')
@click.option('--summary-only', is_flag=True, default=False,
              help='Rebuild the summary TSV from stored results and exit')
@click.option('--debug', is_flag=True, default=False)
def main(db_url, schema, out_dir, base_url, api_key_file, genes, workers,
         max_intron_offset, send_utr_substitutions, batch_size, retry_failed,
         vrs_digest, limit, overwrite, summary_only, debug):
    logging.basicConfig(level=logging.DEBUG if debug else logging.INFO,
                        format='%(asctime)s %(levelname)s %(message)s')
    if max_intron_offset is not None and max_intron_offset < 0:
        max_intron_offset = None    # -1 means "send them anyway"
    if batch_size > DEFAULT_BATCH_SIZE:
        log.warning('Batch size %d is above the %d the ARIANE developers recommend.',
                    batch_size, DEFAULT_BATCH_SIZE)

    results_dir = os.path.join(out_dir, 'results')
    os.makedirs(results_dir, exist_ok=True)
    errors_path = os.path.join(out_dir, 'errors.tsv')
    summary_path = os.path.join(out_dir, 'ariane_summary.tsv')

    if summary_only:
        write_summary(results_dir, summary_path)
        return

    # Before touching the database, so a missing key fails fast.
    api_key = load_api_key(api_key_file)

    gene_list = [g.strip() for g in genes.split(',') if g.strip()]
    conn = psycopg2.connect(db_url, options=f'-c search_path={schema}')
    try:
        query, params = _QUERY, {'genes': gene_list}
        if vrs_digest:
            query += ' AND "VRS_Digest" = %(digest)s'
            params['digest'] = vrs_digest
        with conn.cursor() as cur:
            cur.execute(query, params)
            rows = cur.fetchall()
    finally:
        conn.close()

    if vrs_digest and not rows:
        print(f'Error: VRS digest {vrs_digest!r} not found in the database.')
        return

    errors_lock = threading.Lock()
    errors_fh = open(errors_path, 'w' if overwrite else 'a')
    counts_lock = threading.Lock()
    counts = {'ok': 0, 'done_units': 0}
    error_counts = {}

    def on_error(digest, gene, c_notation, status, message, code=None, retryable=None):
        category = categorize_error(status, message, code, retryable)
        flag = '' if retryable is None else str(bool(retryable)).lower()
        with errors_lock:
            errors_fh.write('\t'.join([
                digest, gene or '', c_notation or '', category, str(status),
                code or '', flag, str(message).replace('\t', ' ')]) + '\n')
            errors_fh.flush()
        with counts_lock:
            error_counts[category] = error_counts.get(category, 0) + 1

    # Normalize notations, and set aside the ones that cannot be sent at all.
    # These go to their own file, rewritten each run: they are a property of the
    # data rather than of this attempt, so appending them to errors.tsv on every
    # resume would just pile up duplicates.
    queryable, unqueryable = [], []
    for digest, gene, hgvs in rows:
        notation = normalize_notation(hgvs)
        if not (gene and notation):
            unqueryable.append((digest, gene, hgvs, 'no gene symbol or no c. notation'))
            continue
        reason = unverifiable_reason(notation, max_intron_offset,
                                     skip_utr_substitutions=not send_utr_substitutions)
        if reason:
            unqueryable.append((digest, gene, notation, reason))
        else:
            queryable.append((digest, gene, notation))
    with open(os.path.join(out_dir, 'unqueryable.tsv'), 'w') as f:
        for digest, gene, notation, reason in unqueryable:
            f.write(f'{digest}\t{gene}\t{notation}\t{reason}\n')
    rows = queryable

    if not overwrite:
        done = load_done(results_dir)
        rows = [r for r in rows if digest_to_filename(r[0])[:-5] not in done]
    if not retry_failed:
        permanent = load_permanent_failures(errors_path)
        if permanent:
            before = len(rows)
            rows = [r for r in rows if r[0] not in permanent]
            print(f'Skipping {before - len(rows)} variant(s) whose earlier failure '
                  f'the service marked non-retryable (--retry-failed resends them).')
    if limit is not None:
        rows = rows[:limit]

    reasons = collections.Counter(u[3].split(',')[0].split(' (')[0] for u in unqueryable)
    for reason, n in reasons.most_common():
        print(f'Not sending {n}: {reason}')

    units = [rows[i:i + batch_size] for i in range(0, len(rows), batch_size)]
    total = len(rows)
    print(f'Classifying {total} variant(s) via {base_url} in batches of {batch_size} '
          f'with {workers} worker(s); {len(unqueryable)} not sent (see unqueryable.tsv).')
    if total == 0:
        write_summary(results_dir, summary_path)
        return

    client = Ariane(base_url, api_key=api_key, debug=debug)
    builds = BuildTracker(api_key=api_key)
    print(f'ARIANE build: {builds.current()[0] or "unknown (could not read /api/resources)"}')
    started = time.time()
    abort = threading.Event()
    auth_failure = {}

    def run_unit(unit):
        if abort.is_set():
            return
        try:
            results = client.classify(unit, on_error)
        except AuthError as e:
            if not abort.is_set():
                auth_failure['message'] = str(e)
                abort.set()
            return
        for digest, gene, c_notation in unit:
            got = results.get(digest)
            if got is None:
                continue
            classification, item_meta = got
            build, checked_at = builds.current(item_meta.get('application_version'))
            payload = result_payload(digest, gene, c_notation, base_url,
                                     classification, item_meta, build, checked_at)
            tmp = os.path.join(results_dir, digest_to_filename(digest) + '.tmp')
            with open(tmp, 'w') as f:
                json.dump(payload, f)
            # Rename so an interrupted write can never leave a half-written file
            # that a resumed run would mistake for a completed variant.
            os.replace(tmp, os.path.join(results_dir, digest_to_filename(digest)))
        with counts_lock:
            counts['ok'] += len(results)
            counts['done_units'] += 1
            if counts['done_units'] % 20 == 0:
                elapsed = time.time() - started
                rate = counts['ok'] / elapsed if elapsed else 0
                print(f'  {counts["ok"]}/{total} classified  ({rate * 3600:.0f}/hour; '
                      f'quota left today: {client.quota_remaining})')
        time.sleep(INTER_REQUEST_DELAY)

    try:
        with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
            list(pool.map(run_unit, units))
    except KeyboardInterrupt:
        print('\nInterrupted; completed variants are on disk and will be skipped '
              'on the next run.')
    finally:
        errors_fh.close()

    elapsed = time.time() - started
    print(f'Done. {counts["ok"]} classified in {elapsed / 3600:.1f}h.')
    if error_counts:
        print('Not classified, by cause:')
        for category, n in sorted(error_counts.items(), key=lambda kv: -kv[1]):
            print(f'  {category:20s} {n}')
    if unqueryable:
        print(f'  {"not sent":20s} {len(unqueryable)}  (see unqueryable.tsv)')
    write_summary(results_dir, summary_path)

    if abort.is_set():
        where = '$ARIANE_API_KEY' if os.environ.get('ARIANE_API_KEY') else api_key_file
        print(f'\nStopped: ARIANE rejected the API key ({auth_failure.get("message")}). '
              f'Check {where}. Variants not yet sent were not recorded as failed, so '
              f'rerunning resumes where this stopped.', file=sys.stderr)
        sys.exit(2)


if __name__ == '__main__':
    main()
