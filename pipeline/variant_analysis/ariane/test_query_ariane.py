"""
Unit tests for query_ariane.py.

The HTTP layer is stubbed: these cover the response handling that decides
whether data is kept or lost -- per-item errors, 422 bisection, 504-as-retry,
auth failure, quota pacing -- without touching the live service or its quota.
"""

import json
import logging
import os
import time

import click
import pytest

import query_ariane as qa


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

V1 = ('ga4gh:VA.one',   'BRCA2', 'c.1975T>G')
V2 = ('ga4gh:VA.two',   'BRCA2', 'c.1978dup')
V3 = ('ga4gh:VA.three', 'BRCA1', 'c.*361G>A')    # the notation ARIANE rejects
V4 = ('ga4gh:VA.four',  'BRCA1', 'c.5219T>G')

REF_ALLELE_MSG = ('Reference allele could not be verified for BRCA1 c.*361 using '
                  'the installed reference-transcript datasets; classification was '
                  'not run.')


class FakeResponse:
    def __init__(self, status_code, body=None, text='', headers=None):
        self.status_code = status_code
        self._body = body
        self.text = text if text else json.dumps(body) if body else ''
        self.headers = headers or {}

    def json(self):
        if self._body is None:
            raise ValueError('no json')
        return self._body


def item_metadata(version='1.9.11'):
    """The per-item metadata v1 returns alongside each classification."""
    return {'api_version': '1.0', 'application_version': version,
            'classifier_fingerprint': 'f' * 64, 'request_id': 'req-1',
            'policy': {'policy_id': 'ENIGMA_BRCA_VCEP_1.2', 'policy_version': '1.2.0'}}


def ok_body(variants):
    """A v1 200 body classifying every variant successfully."""
    return {
        'metadata': {'api_version': '1.0', 'application_version': '1.9.11'},
        'total': len(variants), 'success_count': len(variants), 'error_count': 0,
        'results': [
            {'index': i, 'status': 'ok', 'variant': f'{g} {c}',
             'metadata': item_metadata(),
             'classification': {'predicted_class': 2, 'predicted_label': 'Likely Benign',
                                'total_points': -4, 'c_notation': c}}
            for i, (_, g, c) in enumerate(variants)
        ],
    }


def item_error(index, code, message, retryable=False):
    """A v1 per-item error, as it arrives inside a 200."""
    return {'index': index, 'status': 'error', 'variant': 'x',
            'error': {'code': code, 'message': message, 'retryable': retryable}}


def unprocessable_body(bad_index, msg='Reference allele could not be verified'):
    return {'detail': [{'type': 'value_error',
                        'loc': ['body', 'variants', bad_index],
                        'msg': msg}]}


def auth_body():
    return {'metadata': {'api_version': '1.0'},
            'error': {'code': 'api_key_required',
                      'message': 'Provide an API key in the X-ARIANE-API-Key header',
                      'retryable': False}}


class StubAriane(qa.Ariane):
    """Ariane with _post replaced by a scripted responder."""

    def __init__(self, responder):
        super().__init__('http://stub', api_key='test-key')
        self.responder = responder
        self.calls = []

    def _post(self, variants):
        self.calls.append(list(variants))
        return self.responder(variants)


@pytest.fixture(autouse=True)
def no_sleep(monkeypatch):
    """Backoff and quota sleeps would make these tests take minutes."""
    monkeypatch.setattr(qa.time, 'sleep', lambda *_: None)


def collect_errors():
    errs = []

    def on_error(d, g, c, s, m, code=None, retryable=None):
        errs.append((d, s, str(m), code, retryable))
    return errs, on_error


# ---------------------------------------------------------------------------
# which variants are sent at all
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('notation,offset', [
    ('c.5074+1G>A', 1),
    ('c.9257-1532_9257-1530del', 1532),
    ('c.682-24G>A', 24),
    ('c.1975T>G', None),        # coding, no intron offset
    ('c.*361G>A', None),        # 3' UTR
    ('', None),
])
def test_intron_offset(notation, offset):
    assert qa.intron_offset(notation) == offset


@pytest.mark.parametrize('notation,reason', [
    # ARIANE verifies the CDS plus 50 bases of intronic flank. Substitutions
    # beyond that are a guaranteed rejection, so they are not sent.
    ('c.5074+1G>A', None),
    ('c.5074+50G>A', None),     # exactly at the boundary -- still classifiable
    ('c.5074+51G>A', 'deep intronic'),   # one past it
    ('c.682-7256G>A', 'deep intronic'),
    # the bundle holds no UTR sequence, so UTR substitutions are rejected too
    ('c.*361G>A', 'UTR substitution'),
    ('c.-20+107G>A', 'UTR substitution'),
    # indels carry no reference allele to verify, so depth does not matter
    ('c.9257-1532_9257-1530del', None),
    ('c.5074+900dup', None),
    ('c.*361dup', None),
    ('c.-57_-55del', None),
    ('c.1975T>G', None),
])
def test_unverifiable_reason(notation, reason):
    got = qa.unverifiable_reason(notation)
    if reason is None:
        assert got is None
    else:
        assert got is not None and got.startswith(reason)


def test_utr_substitutions_can_be_sent_on_request():
    """--send-utr-substitutions turns off only the UTR rule."""
    assert qa.unverifiable_reason('c.*361G>A', skip_utr_substitutions=False) is None
    assert qa.unverifiable_reason('c.682-7256G>A', skip_utr_substitutions=False)


def test_unverifiable_reason_disabled_sends_deep_intronic():
    """--max-intron-offset -1 turns off the intronic rule."""
    assert qa.unverifiable_reason('c.682-7256G>A', max_intron_offset=None) is None


def test_digest_to_filename_sanitizes_colon():
    assert qa.digest_to_filename('ga4gh:VA.abc') == 'ga4gh_VA.abc.json'


@pytest.mark.parametrize('status,message,code,retryable,expected', [
    # pre-v1 shapes, still recognised
    (422, REF_ALLELE_MSG, None, None, 'service_limitation'),
    (422, 'some other validation problem', None, None, 'error'),
    ('skipped', 'no gene symbol or no c. notation', None, None, 'data'),
    ('missing', 'no result at this index', None, None, 'data'),
    ('retries_exhausted', 'no success after 4 attempts', None, None, 'transient'),
    (500, 'internal server error', None, None, 'error'),
    ('error', None, None, None, 'error'),
    # v1 per-item errors
    ('error', REF_ALLELE_MSG, 'invalid_variant', False, 'service_limitation'),
    ('error', 'Invalid c. HGVS for BRCA1: c.NOT_HGVS', 'invalid_variant', False, 'data'),
    ('error', 'no population evidence', 'structural_population_evidence_unavailable',
     False, 'service_limitation'),
    ('error', 'SpliceAI is required ... (no_grch38_coords)',
     'spliceai_coordinates_unavailable', False, 'service_limitation'),
    ('error', 'SpliceAI is required for this automatic classification but is '
              'unavailable (api_error): SpliceAI lookup timed out after 30 seconds.',
     'spliceai_unavailable', False, 'transient'),
    ('error', 'upstream busy', 'busy', True, 'transient'),
])
def test_categorize_error(status, message, code, retryable, expected):
    assert qa.categorize_error(status, message, code, retryable) == expected


@pytest.mark.parametrize('stored,expected', [
    ('c.1975T>G', 'c.1975T>G'),
    # transcript-qualified rows are well-formed, just prefixed
    ('NM_007294.3:c.*6207C>T', 'c.*6207C>T'),
    ('NM_007294.3:c.5219T>G', 'c.5219T>G'),
    ('  c.68_69delAG  ', 'c.68_69delAG'),
    ('-', None),               # the table's placeholder for missing
    ('', None),
    (None, None),
    ('p.Val1740Gly', None),    # protein notation is not a c. notation
])
def test_normalize_notation(stored, expected):
    assert qa.normalize_notation(stored) == expected


# ---------------------------------------------------------------------------
# the API key
# ---------------------------------------------------------------------------

def key_file(tmp_path, content='ariane_v1_test\n', mode=0o600):
    path = tmp_path / 'api_key'
    path.write_text(content)
    os.chmod(path, mode)
    return str(path)


def test_api_key_env_takes_precedence(monkeypatch, tmp_path):
    monkeypatch.setenv('ARIANE_API_KEY', 'from-env')
    assert qa.load_api_key(key_file(tmp_path)) == 'from-env'


def test_api_key_file_is_read_and_stripped(monkeypatch, tmp_path):
    monkeypatch.delenv('ARIANE_API_KEY', raising=False)
    assert qa.load_api_key(key_file(tmp_path)) == 'ariane_v1_test'


def test_missing_api_key_names_the_file(monkeypatch, tmp_path):
    monkeypatch.delenv('ARIANE_API_KEY', raising=False)
    path = str(tmp_path / 'nope')
    with pytest.raises(click.ClickException) as exc:
        qa.load_api_key(path)
    assert path in exc.value.message


def test_empty_api_key_file_is_an_error(monkeypatch, tmp_path):
    monkeypatch.delenv('ARIANE_API_KEY', raising=False)
    with pytest.raises(click.ClickException):
        qa.load_api_key(key_file(tmp_path, content='\n'))


def test_readable_key_file_warns(monkeypatch, tmp_path, caplog):
    monkeypatch.delenv('ARIANE_API_KEY', raising=False)
    with caplog.at_level(logging.WARNING):
        qa.load_api_key(key_file(tmp_path, mode=0o644))
    assert 'chmod 600' in caplog.text
    assert 'ariane_v1_test' not in caplog.text, 'the key must never be logged'


def test_key_is_sent_as_the_ariane_header():
    client = qa.Ariane('http://stub', api_key='k-123')
    assert client.session.headers[qa.API_KEY_HEADER] == 'k-123'


def test_no_key_header_without_a_key():
    assert qa.API_KEY_HEADER not in qa.Ariane('http://stub').session.headers


# ---------------------------------------------------------------------------
# happy path and per-item errors
# ---------------------------------------------------------------------------

def test_successful_batch_returns_classification_and_metadata():
    variants = [V1, V2, V4]
    client = StubAriane(lambda v: FakeResponse(200, ok_body(v)))
    errs, on_error = collect_errors()
    out = client.classify(variants, on_error)
    assert set(out) == {V1[0], V2[0], V4[0]}
    classification, metadata = out[V1[0]]
    assert classification['predicted_label'] == 'Likely Benign'
    assert metadata['application_version'] == '1.9.11'
    assert errs == []


def test_pre_v1_result_key_is_still_understood():
    body = ok_body([V1])
    body['results'][0]['result'] = body['results'][0].pop('classification')
    client = StubAriane(lambda v: FakeResponse(200, body))
    errs, on_error = collect_errors()
    out = client.classify([V1], on_error)
    assert out[V1[0]][0]['predicted_label'] == 'Likely Benign'


def test_per_item_error_is_recorded_with_code_and_retryable():
    """v1 reports a failed item inside a 200; its siblings must survive, and the
    error's code and retryable flag must reach the error log."""
    variants = [V1, V2]
    body = ok_body(variants)
    body['results'][1] = item_error(1, 'invalid_variant', REF_ALLELE_MSG, False)
    client = StubAriane(lambda v: FakeResponse(200, body))
    errs, on_error = collect_errors()
    out = client.classify(variants, on_error)
    assert set(out) == {V1[0]}
    assert errs == [(V2[0], 'error', REF_ALLELE_MSG, 'invalid_variant', False)]


def test_legacy_string_error_is_still_recorded():
    variants = [V1, V2]
    body = ok_body(variants)
    body['results'][1] = {'index': 1, 'status': 'error', 'variant': 'x',
                          'error': 'could not normalise'}
    client = StubAriane(lambda v: FakeResponse(200, body))
    errs, on_error = collect_errors()
    client.classify(variants, on_error)
    assert errs[0][0] == V2[0] and 'could not normalise' in errs[0][2]


def test_missing_index_in_response_is_an_error_not_a_silent_drop():
    variants = [V1, V2]
    body = ok_body(variants)
    del body['results'][1]
    client = StubAriane(lambda v: FakeResponse(200, body))
    errs, on_error = collect_errors()
    out = client.classify(variants, on_error)
    assert set(out) == {V1[0]}
    assert errs[0][0] == V2[0] and errs[0][1] == 'missing'


# ---------------------------------------------------------------------------
# a rejected key stops the run instead of draining the queue
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('status', [401, 403])
def test_rejected_key_raises_and_records_nothing(status):
    """Last time a missing key cost 7,016 requests, each recorded as a failure.
    Now the first rejection stops the batch, and no variant is marked failed,
    so a rerun with a working key picks them all up."""
    client = StubAriane(lambda v: FakeResponse(status, auth_body()))
    errs, on_error = collect_errors()
    with pytest.raises(qa.AuthError) as exc:
        client.classify([V1, V2], on_error)
    assert 'api_key_required' in str(exc.value)
    assert errs == []
    assert len(client.calls) == 1


# ---------------------------------------------------------------------------
# 422: one bad variant must not take the batch down with it
# ---------------------------------------------------------------------------

def test_422_bisects_and_keeps_the_healthy_variants():
    """v1 reports bad variants per item, but should a whole payload ever be
    rejected again, bisection must isolate the offender and keep the rest."""
    variants = [V1, V2, V3, V4]      # V3 is the bad one

    def responder(v):
        if any(x[0] == V3[0] for x in v):
            bad = [i for i, x in enumerate(v) if x[0] == V3[0]][0]
            return FakeResponse(422, unprocessable_body(bad))
        return FakeResponse(200, ok_body(v))

    client = StubAriane(responder)
    errs, on_error = collect_errors()
    out = client.classify(variants, on_error)

    assert set(out) == {V1[0], V2[0], V4[0]}, 'healthy variants must survive'
    assert [e[0] for e in errs] == [V3[0]], 'only the offender is an error'
    assert errs[0][1] == 422
    assert 'Reference allele could not be verified' in errs[0][2]


def test_422_on_a_singleton_records_the_variant_and_stops():
    client = StubAriane(lambda v: FakeResponse(422, unprocessable_body(0)))
    errs, on_error = collect_errors()
    out = client.classify([V3], on_error)
    assert out == {}
    assert [e[0] for e in errs] == [V3[0]]
    # A singleton cannot be split further, so exactly one request is made.
    assert len(client.calls) == 1


def test_422_with_several_bad_variants_isolates_all_of_them():
    variants = [V1, V3, V4, ('ga4gh:VA.five', 'BRCA1', 'c.*99A>T')]
    bad = {V3[0], 'ga4gh:VA.five'}

    def responder(v):
        idx = [i for i, x in enumerate(v) if x[0] in bad]
        if idx:
            return FakeResponse(422, unprocessable_body(idx[0]))
        return FakeResponse(200, ok_body(v))

    client = StubAriane(responder)
    errs, on_error = collect_errors()
    out = client.classify(variants, on_error)
    assert set(out) == {V1[0], V4[0]}
    assert {e[0] for e in errs} == bad


# ---------------------------------------------------------------------------
# 504: the backend finished anyway, so retry rather than fail
# ---------------------------------------------------------------------------

def test_504_is_retried_and_the_cached_result_collected():
    variants = [V1, V2]
    state = {'n': 0}

    def responder(v):
        state['n'] += 1
        if state['n'] == 1:
            return FakeResponse(504, text='gateway timeout')
        return FakeResponse(200, ok_body(v))

    client = StubAriane(responder)
    errs, on_error = collect_errors()
    out = client.classify(variants, on_error)
    assert set(out) == {V1[0], V2[0]}, 'a 504 must not lose the batch'
    assert errs == []
    assert len(client.calls) == 2


def test_persistent_failure_eventually_reports_every_variant():
    client = StubAriane(lambda v: FakeResponse(504, text='gateway timeout'))
    errs, on_error = collect_errors()
    out = client.classify([V1, V2], on_error)
    assert out == {}
    assert {e[0] for e in errs} == {V1[0], V2[0]}
    assert all(e[1] == 'retries_exhausted' for e in errs)
    assert len(client.calls) == qa.MAX_ATTEMPTS


def test_connection_error_is_retried():
    import requests as rq
    state = {'n': 0}

    def responder(v):
        state['n'] += 1
        if state['n'] == 1:
            raise rq.ConnectionError('connection reset')
        return FakeResponse(200, ok_body(v))

    client = StubAriane(responder)
    errs, on_error = collect_errors()
    out = client.classify([V1], on_error)
    assert set(out) == {V1[0]}
    assert errs == []


def test_unexpected_status_is_not_retried():
    client = StubAriane(lambda v: FakeResponse(404, text='not found'))
    errs, on_error = collect_errors()
    out = client.classify([V1], on_error)
    assert out == {}
    assert errs[0][1] == 404
    assert len(client.calls) == 1, 'a 404 will not fix itself; do not hammer'


# ---------------------------------------------------------------------------
# quota: 5,000 classifications per key per UTC day
# ---------------------------------------------------------------------------

def test_quota_headers_are_recorded():
    headers = {'x-ratelimit-remaining': '42', 'x-ratelimit-reset': '1789171200'}
    client = StubAriane(lambda v: FakeResponse(200, ok_body(v), headers=headers))
    client.classify([V1], collect_errors()[1])
    assert client.quota_remaining == 42
    assert client.quota_reset == 1789171200


def test_short_quota_sleeps_until_the_reset(monkeypatch):
    slept = []
    monkeypatch.setattr(qa.time, 'sleep', slept.append)
    client = qa.Ariane('http://stub')
    client.quota_remaining, client.quota_reset = 2, int(time.time()) + 3600
    client.wait_for_quota(5)
    assert len(slept) == 1
    assert 3600 <= slept[0] <= 3600 + qa.QUOTA_RESET_MARGIN + 5
    assert client.quota_remaining is None, 'a new window; the next response says'


def test_enough_quota_does_not_sleep(monkeypatch):
    slept = []
    monkeypatch.setattr(qa.time, 'sleep', slept.append)
    client = qa.Ariane('http://stub')
    client.quota_remaining, client.quota_reset = 5, int(time.time()) + 3600
    client.wait_for_quota(5)
    assert slept == []


def test_quota_429_waits_for_the_reset_without_spending_attempts(monkeypatch):
    """Exhausting the daily quota is not a failure. More 429s than MAX_ATTEMPTS
    must still end in success once the window resets."""
    slept = []
    monkeypatch.setattr(qa.time, 'sleep', slept.append)
    reset = int(time.time()) + 600
    spent = {'x-ratelimit-remaining': '0', 'x-ratelimit-reset': str(reset)}
    state = {'n': 0}

    def responder(v):
        state['n'] += 1
        if state['n'] <= qa.MAX_ATTEMPTS + 2:
            return FakeResponse(429, text='quota', headers=spent)
        return FakeResponse(200, ok_body(v))

    client = StubAriane(responder)
    errs, on_error = collect_errors()
    out = client.classify([V1], on_error)
    assert set(out) == {V1[0]}
    assert errs == []
    assert max(slept) >= 600 - 5, 'it waited for the reset, not a short backoff'


def test_rate_limit_429_uses_retry_after():
    state = {'n': 0}

    def responder(v):
        state['n'] += 1
        if state['n'] == 1:
            return FakeResponse(429, text='slow down', headers={'Retry-After': '1'})
        return FakeResponse(200, ok_body(v))

    client = StubAriane(responder)
    errs, on_error = collect_errors()
    out = client.classify([V1], on_error)
    assert set(out) == {V1[0]}
    assert errs == []


def test_endless_429s_eventually_give_up(monkeypatch):
    monkeypatch.setattr(qa, 'MAX_CONSECUTIVE_429', 3)
    client = StubAriane(lambda v: FakeResponse(429, text='slow down',
                                               headers={'Retry-After': '1'}))
    errs, on_error = collect_errors()
    assert client.classify([V1], on_error) == {}
    assert errs[0][1] == 'retries_exhausted'
    assert len(client.calls) == 3


# ---------------------------------------------------------------------------
# build version
# ---------------------------------------------------------------------------

def test_build_is_refreshed_on_version_change_and_hourly(monkeypatch):
    tracker = qa.BuildTracker()
    fetches = []
    monkeypatch.setattr(tracker, '_fetch', lambda: fetches.append(1) or '626d234')
    clock = [1000.0]
    monkeypatch.setattr(qa.time, 'monotonic', lambda: clock[0])

    assert tracker.current()[0] == '626d234'
    assert len(fetches) == 1
    tracker.current('1.9.11')            # first version seen: nothing to compare
    tracker.current('1.9.11')
    assert len(fetches) == 1, 'no refetch while nothing changes'
    tracker.current('1.9.12')            # redeployed
    assert len(fetches) == 2
    clock[0] += qa.BUILD_REFRESH_SECONDS
    tracker.current('1.9.12')            # an hour on
    assert len(fetches) == 3


def test_build_keeps_last_known_value_when_a_refresh_fails(monkeypatch):
    tracker = qa.BuildTracker()
    answers = iter(['626d234', None])
    monkeypatch.setattr(tracker, '_fetch', lambda: next(answers))
    tracker.current('1.9.11')
    build, checked_at = tracker.current('1.9.12')    # refresh fails
    assert build == '626d234'
    assert checked_at is not None


def test_build_fetch_sends_the_api_key(monkeypatch):
    """/api/resources began returning 401 without a key on 2026-09-11; without
    the header every result would carry build_version None."""
    seen = {}

    class Resp:
        def raise_for_status(self):
            pass

        def json(self):
            return {'build_revision': '626d234'}

    def fake_get(url, headers=None, timeout=None):
        seen['headers'] = headers
        return Resp()

    monkeypatch.setattr(qa.requests, 'get', fake_get)
    assert qa.BuildTracker(api_key='k-123')._fetch() == '626d234'
    assert seen['headers'] == {qa.API_KEY_HEADER: 'k-123'}


def test_result_payload_stamps_build_and_keeps_ariane_metadata():
    payload = qa.result_payload(
        'ga4gh:VA.one', 'BRCA2', 'c.1975T>G', 'http://stub/api/v1',
        {'predicted_class': 2}, item_metadata(), '626d234', '2026-09-11T18:00:00+00:00')
    ariane = payload['_meta']['ariane']
    assert ariane['application_version'] == '1.9.11'
    assert ariane['classifier_fingerprint'] == 'f' * 64
    assert ariane['build_version'] == '626d234'
    assert ariane['build_version_checked_at'] == '2026-09-11T18:00:00+00:00'
    assert payload['result'] == {'predicted_class': 2}


# ---------------------------------------------------------------------------
# resumability
# ---------------------------------------------------------------------------

def test_load_done_reads_completed_digests(tmp_path):
    results = tmp_path / 'results'
    results.mkdir()
    (results / 'ga4gh_VA.one.json').write_text('{}')
    (results / 'ga4gh_VA.two.json').write_text('{}')
    (results / 'ga4gh_VA.three.json.tmp').write_text('{}')   # partial write
    done = qa.load_done(str(results))
    assert done == {'ga4gh_VA.one', 'ga4gh_VA.two'}, \
        'a .tmp file is an interrupted write, not a completed variant'


def test_permanent_failures_are_skipped_on_resume(tmp_path):
    """Under a daily quota, resending known-permanent failures on every restart
    spends real classifications; transient ones are worth another try."""
    errors = tmp_path / 'errors.tsv'
    errors.write_text('\n'.join([
        # non-retryable service limitation: skip
        'ga4gh:VA.a\tBRCA1\tc.*1A>G\tservice_limitation\terror\tinvalid_variant\tfalse\t' + REF_ALLELE_MSG,
        # SpliceAI timeout, flagged non-retryable but transient in fact: retry
        'ga4gh:VA.b\tBRCA2\tc.1A>G\ttransient\terror\tspliceai_unavailable\tfalse\t(api_error) timed out',
        # retryable flag unknown: retry
        'ga4gh:VA.c\tBRCA2\tc.2A>G\terror\t500\t\t\tboom',
        # a row from before errors.tsv carried the retryable column: ignore
        'ga4gh:VA.d\tBRCA2\tc.3A>G\tservice_limitation\t422\told message',
    ]) + '\n')
    assert qa.load_permanent_failures(str(errors)) == {'ga4gh:VA.a'}
    assert qa.load_permanent_failures(str(tmp_path / 'absent.tsv')) == set()


def test_write_summary_reads_stored_results(tmp_path):
    results = tmp_path / 'results'
    results.mkdir()
    (results / 'ga4gh_VA.one.json').write_text(json.dumps({
        '_meta': {'VRS_Digest': 'ga4gh:VA.one', 'gene': 'BRCA2',
                  'c_notation': 'c.1975T>G',
                  'ariane': {'application_version': '1.9.11',
                             'build_version': '626d234'}},
        # criteria come back as {name, applies, strength, points, reason, ...};
        # only the ones that applied belong in the summary.
        'result': {'predicted_class': 2, 'predicted_label': 'Likely Benign',
                   'total_points': -4, 'evidence_direction': 'benign',
                   'criteria': [
                       {'name': 'BS1', 'applies': True, 'strength': 'Strong'},
                       {'name': 'BP4', 'applies': True, 'strength': 'Supporting'},
                       {'name': 'PM2', 'applies': False, 'strength': 'Moderate'},
                   ]},
    }))
    summary = tmp_path / 'summary.tsv'
    qa.write_summary(str(results), str(summary))
    lines = summary.read_text().strip().split('\n')
    assert lines[0].startswith('VRS_Digest\tgene\tc_notation\tpredicted_class')
    assert lines[0].endswith('ariane_version\tbuild_version')
    assert 'ga4gh:VA.one' in lines[1]
    assert 'Likely Benign' in lines[1]
    assert 'BS1(Strong);BP4(Supporting)' in lines[1]
    assert 'PM2' not in lines[1], 'criteria that did not apply are not "met"'
    assert lines[1].endswith('1.9.11\t626d234')
