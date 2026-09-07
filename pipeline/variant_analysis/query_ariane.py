"""
Classify every variant with the ARIANE service, storing results on the filesystem.

ARIANE (https://ariane-app.duckdns.org/) is an ENIGMA BRCA1/2 ACMG
classification service. This script reads gene + HGVS cDNA from the variant
table, POSTs them to /api/classify/batch, and writes one JSON file per variant.
Nothing is written to the database.

Output layout (under --out-dir):
    results/<VRS_Digest with ':' replaced by '_'>.json
    errors.tsv                variants the service could not classify
    ariane_summary.tsv        one row per classified variant
    run.log                   (when launched detached)

Variants that already have a result file are skipped unless --overwrite, so the
run is idempotent and a multi-day pass can be stopped and resumed freely.

Three measured properties of the service shape this script:

  * A batch is processed serially server-side at ~5.3 s per uncached variant,
    and nginx cuts requests off at 60 s. Batches of 5 (~26 s) leave margin;
    a batch of 20 reliably 504s.
  * A 504 does not mean the work was lost. The backend completes and caches it,
    so re-POSTing the same batch returns quickly. 504 is retried, not failed.
  * A single unclassifiable variant rejects the WHOLE batch with 422. Variants
    whose notation is likely to be rejected are therefore sent individually, and
    any batch that still 422s is bisected to isolate the offender rather than
    losing its healthy siblings.
"""

import concurrent.futures
import datetime
import json
import logging
import os
import random
import threading
import time

import click
import psycopg2
import requests

log = logging.getLogger(__name__)

DEFAULT_BASE_URL = 'https://ariane-app.duckdns.org'

# Batch of 5 takes ~26s server-side against a 60s gateway timeout.
DEFAULT_BATCH_SIZE = 5
DEFAULT_WORKERS = 2

REQUEST_TIMEOUT = 90        # > the 60s gateway timeout, so we see the 504 itself
MAX_ATTEMPTS = 4
INTER_REQUEST_DELAY = 0.5   # politeness: this is a small shared instance

_QUERY = """
SELECT "VRS_Digest", "Gene_Symbol", "HGVS_cDNA"
FROM variant
WHERE "Gene_Symbol" = ANY(%(genes)s)
"""

def categorize_error(status, message):
    """Bucket a failure so the ~2,800 expected ones do not hide real breakage.

    ARIANE's installed reference bundle carries CDS sequence only, so any UTR
    *substitution* fails reference-allele verification (UTR indels pass, since
    they need no reference base). Those are a known limit of the service, not a
    fault of this run, and are labelled as such.
    """
    text = str(message or '')
    if status == 422 and 'Reference allele could not be verified' in text:
        return 'service_limitation'
    if status in ('skipped', 'missing'):
        return 'data'
    if status == 'retries_exhausted':
        return 'transient'
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


# Notations in these shapes are the ones ARIANE is most likely to reject (the
# confirmed case was a 3' UTR variant whose reference allele could not be
# verified). They are sent one at a time so a rejection cannot take a batch of
# healthy variants down with it.
def is_risky_notation(c_notation):
    return not (c_notation or '').startswith('c.') or \
           c_notation.startswith('c.*') or c_notation.startswith('c.-')


def digest_to_filename(vrs_digest):
    """Filesystem-safe filename for a VRS digest."""
    return vrs_digest.replace(':', '_') + '.json'


class Ariane:
    """Client for the ARIANE batch classification endpoint."""

    def __init__(self, base_url, debug=False):
        self.base_url = base_url.rstrip('/')
        self.debug = debug
        self._local = threading.local()

    @property
    def session(self):
        # One session per worker thread; requests.Session is not thread-safe.
        if not hasattr(self._local, 'session'):
            self._local.session = requests.Session()
        return self._local.session

    def _post(self, variants):
        payload = {'variants': [{'gene': g, 'c_notation': c} for _, g, c in variants]}
        return self.session.post(f'{self.base_url}/api/classify/batch',
                                 json=payload, timeout=REQUEST_TIMEOUT)

    def classify(self, variants, on_error):
        """Classify a list of (digest, gene, c_notation).

        Returns {digest: result_dict}. Variants the service rejects are passed
        to on_error(digest, gene, c_notation, status, message) instead.
        """
        results = {}
        for attempt in range(1, MAX_ATTEMPTS + 1):
            try:
                resp = self._post(variants)
            except requests.RequestException as e:
                self._backoff(attempt, f'{type(e).__name__}: {e}')
                continue

            if resp.status_code == 200:
                return self._collect(variants, resp.json(), on_error)

            if resp.status_code == 422:
                # The whole payload was rejected because of at least one bad
                # variant. Split to find it; the healthy ones are cached by now,
                # so the extra round trips are cheap.
                return self._bisect(variants, resp, on_error)

            if resp.status_code == 504:
                # The backend finished and cached the work even though nginx
                # gave up; retrying returns it quickly.
                log.info('504 on %d variant(s); retrying to collect cached results',
                         len(variants))
                self._backoff(attempt, '504 gateway timeout', base=2.0)
                continue

            if resp.status_code in (429, 500, 502, 503):
                retry_after = resp.headers.get('Retry-After')
                delay = float(retry_after) if (retry_after or '').isdigit() else None
                self._backoff(attempt, f'HTTP {resp.status_code}', fixed=delay)
                continue

            # Anything else is not worth retrying.
            for digest, gene, c in variants:
                on_error(digest, gene, c, resp.status_code, resp.text[:300])
            return results

        for digest, gene, c in variants:
            on_error(digest, gene, c, 'retries_exhausted',
                     f'no success after {MAX_ATTEMPTS} attempts')
        return results

    def _collect(self, variants, body, on_error):
        """Pull per-item results out of a 200 response."""
        results = {}
        by_index = {r.get('index'): r for r in body.get('results', [])}
        for i, (digest, gene, c) in enumerate(variants):
            item = by_index.get(i)
            if item is None:
                on_error(digest, gene, c, 'missing', 'no result at this index')
            elif item.get('status') == 'ok' and item.get('result') is not None:
                results[digest] = item['result']
            else:
                on_error(digest, gene, c, item.get('status', 'error'),
                         item.get('error') or 'no result returned')
        return results

    def _bisect(self, variants, resp, on_error):
        """A 422 rejects the entire payload. Split until the offender is alone."""
        if len(variants) == 1:
            digest, gene, c = variants[0]
            on_error(digest, gene, c, 422, self._422_message(resp))
            return {}
        mid = len(variants) // 2
        out = {}
        for half in (variants[:mid], variants[mid:]):
            out.update(self.classify(half, on_error))
        return out

    @staticmethod
    def _422_message(resp):
        try:
            detail = resp.json().get('detail')
            if isinstance(detail, list) and detail:
                return detail[0].get('msg', '')[:300]
            return str(detail)[:300]
        except ValueError:
            return resp.text[:300]

    def _backoff(self, attempt, reason, base=1.5, fixed=None):
        delay = fixed if fixed is not None else base * (2 ** (attempt - 1))
        delay += random.uniform(0, 0.5)     # jitter, so workers do not sync up
        log.warning('attempt %d failed (%s); sleeping %.1fs', attempt, reason, delay)
        time.sleep(delay)


def load_done(results_dir):
    """Digests already fetched — the run's checkpoint."""
    return {f[:-len('.json')] for f in os.listdir(results_dir) if f.endswith('.json')}


def write_summary(results_dir, summary_path):
    """One row per classified variant, from the stored JSON."""
    import csv
    rows = []
    for name in sorted(os.listdir(results_dir)):
        if not name.endswith('.json'):
            continue
        with open(os.path.join(results_dir, name)) as f:
            d = json.load(f)
        meta, res = d.get('_meta', {}), d.get('result', {})
        # Each criterion is {name, applies, strength, points, reason, ...};
        # only the ones that actually applied are worth carrying here.
        applied = [c for c in (res.get('criteria') or []) if c.get('applies')]
        rows.append([
            meta.get('VRS_Digest'), meta.get('gene'), meta.get('c_notation'),
            res.get('predicted_class'), res.get('predicted_label'),
            res.get('total_points'), res.get('evidence_direction'),
            ';'.join(f'{c.get("name")}({c.get("strength")})' for c in applied),
        ])
    with open(summary_path, 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['VRS_Digest', 'gene', 'c_notation', 'predicted_class',
                    'predicted_label', 'total_points', 'evidence_direction',
                    'criteria_met'])
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
@click.option('--genes', default='BRCA1,BRCA2', show_default=True)
@click.option('--workers', default=DEFAULT_WORKERS, show_default=True,
              help='Concurrent requests. Kept low deliberately: this is a small '
                   'shared instance, and the run is resumable so speed is cheap '
                   'to trade away.')
@click.option('--batch-size', default=DEFAULT_BATCH_SIZE, show_default=True,
              help='Variants per request. Above ~8 the 60s gateway timeout hits.')
@click.option('--vrs-digest', default=None, metavar='DIGEST',
              help='Classify only this variant')
@click.option('--limit', default=None, type=int, help='Classify at most this many')
@click.option('--overwrite', is_flag=True, default=False,
              help='Re-fetch variants that already have a result file')
@click.option('--summary-only', is_flag=True, default=False,
              help='Rebuild the summary TSV from stored results and exit')
@click.option('--debug', is_flag=True, default=False)
def main(db_url, schema, out_dir, base_url, genes, workers, batch_size,
         vrs_digest, limit, overwrite, summary_only, debug):
    logging.basicConfig(level=logging.DEBUG if debug else logging.INFO,
                        format='%(asctime)s %(levelname)s %(message)s')

    results_dir = os.path.join(out_dir, 'results')
    os.makedirs(results_dir, exist_ok=True)
    errors_path = os.path.join(out_dir, 'errors.tsv')
    summary_path = os.path.join(out_dir, 'ariane_summary.tsv')

    if summary_only:
        write_summary(results_dir, summary_path)
        return

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

    def on_error(digest, gene, c_notation, status, message):
        category = categorize_error(status, message)
        with errors_lock:
            errors_fh.write(f'{digest}\t{gene}\t{c_notation}\t{category}\t{status}\t'
                            f'{str(message).replace(chr(9), " ")}\n')
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
        if gene and notation:
            queryable.append((digest, gene, notation))
        else:
            unqueryable.append((digest, gene, hgvs))
    with open(os.path.join(out_dir, 'unqueryable.tsv'), 'w') as f:
        for digest, gene, hgvs in unqueryable:
            f.write(f'{digest}\t{gene}\t{hgvs}\tno gene symbol or no c. notation\n')
    rows = queryable

    if not overwrite:
        done = load_done(results_dir)
        rows = [r for r in rows if digest_to_filename(r[0])[:-5] not in done]
    if limit is not None:
        rows = rows[:limit]

    # Risky notations go one at a time so a 422 cannot take healthy variants
    # down with it; the rest are batched.
    risky = [r for r in rows if is_risky_notation(r[2])]
    safe = [r for r in rows if not is_risky_notation(r[2])]
    units = ([[r] for r in risky]
             + [safe[i:i + batch_size] for i in range(0, len(safe), batch_size)])
    random.shuffle(units)   # spread risky singletons through the run

    total = len(rows)
    print(f'Classifying {total} variant(s) via {base_url} '
          f'({len(safe)} batched {batch_size}/request, {len(risky)} sent singly) '
          f'with {workers} worker(s); {len(unqueryable)} unqueryable skipped.')
    if total == 0:
        write_summary(results_dir, summary_path)
        return

    client = Ariane(base_url, debug=debug)
    started = time.time()

    def run_unit(unit):
        results = client.classify(unit, on_error)
        for digest, gene, c_notation in unit:
            res = results.get(digest)
            if res is None:
                continue
            payload = {
                '_meta': {
                    'VRS_Digest': digest, 'gene': gene, 'c_notation': c_notation,
                    'source': base_url,
                    'retrieved': datetime.datetime.now(
                        datetime.timezone.utc).isoformat(timespec='seconds'),
                },
                'result': res,
            }
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
                remaining = (total - counts['ok']) / rate / 3600 if rate else 0
                print(f'  {counts["ok"]}/{total} classified  '
                      f'({rate * 3600:.0f}/hour, ~{remaining:.1f}h remaining)')
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
            note = ('  (ARIANE reference bundle is CDS-only, so UTR '
                    'substitutions cannot be verified)'
                    if category == 'service_limitation' else '')
            print(f'  {category:20s} {n}{note}')
    if unqueryable:
        print(f'  {"unqueryable":20s} {len(unqueryable)}'
              '  (no gene or no c. notation; see unqueryable.tsv)')
    write_summary(results_dir, summary_path)


if __name__ == '__main__':
    main()
