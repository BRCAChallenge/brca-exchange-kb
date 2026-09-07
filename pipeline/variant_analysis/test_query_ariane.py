"""
Unit tests for query_ariane.py.

The HTTP layer is stubbed: these cover the response handling that decides
whether data is kept or lost -- 422 bisection, 504-as-retry, per-item errors --
without touching the live service.
"""

import json

import pytest

import query_ariane as qa


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

V1 = ('ga4gh:VA.one',   'BRCA2', 'c.1975T>G')
V2 = ('ga4gh:VA.two',   'BRCA2', 'c.1978dup')
V3 = ('ga4gh:VA.three', 'BRCA1', 'c.*361G>A')    # the notation ARIANE rejects
V4 = ('ga4gh:VA.four',  'BRCA1', 'c.5219T>G')


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


def ok_body(variants):
    """A 200 body classifying every variant successfully."""
    return {
        'total': len(variants), 'success_count': len(variants), 'error_count': 0,
        'results': [
            {'index': i, 'status': 'ok', 'variant': f'{g} {c}',
             'error': None,
             'result': {'predicted_class': 2, 'predicted_label': 'Likely Benign',
                        'total_points': -4, 'c_notation': c}}
            for i, (_, g, c) in enumerate(variants)
        ],
    }


def unprocessable_body(bad_index, msg='Reference allele could not be verified'):
    return {'detail': [{'type': 'value_error',
                        'loc': ['body', 'variants', bad_index],
                        'msg': msg}]}


class StubAriane(qa.Ariane):
    """Ariane with _post replaced by a scripted responder."""

    def __init__(self, responder):
        super().__init__('http://stub')
        self.responder = responder
        self.calls = []

    def _post(self, variants):
        self.calls.append(list(variants))
        return self.responder(variants)


@pytest.fixture(autouse=True)
def no_sleep(monkeypatch):
    """Backoff sleeps would make these tests take minutes."""
    monkeypatch.setattr(qa.time, 'sleep', lambda *_: None)


def collect_errors():
    errs = []
    return errs, lambda d, g, c, s, m: errs.append((d, s, str(m)))


# ---------------------------------------------------------------------------
# notation routing
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('notation,risky', [
    ('c.1975T>G', False),
    ('c.68_69delAG', False),
    ('c.5074+1G>A', False),      # intronic is fine
    ('c.*361G>A', True),         # 3' UTR -- the confirmed rejection
    ('c.-20+107dup', True),      # 5' UTR
    ('n.1468A>T', True),         # not a c. notation at all
    ('', True),
    (None, True),
])
def test_is_risky_notation(notation, risky):
    assert qa.is_risky_notation(notation) is risky


def test_digest_to_filename_sanitizes_colon():
    assert qa.digest_to_filename('ga4gh:VA.abc') == 'ga4gh_VA.abc.json'


@pytest.mark.parametrize('status,message,expected', [
    # The ~2,800 expected UTR-substitution failures must be distinguishable
    # from real breakage in the error log.
    (422, 'Reference allele could not be verified for BRCA1 c.*12 using the '
          'installed reference-transcript datasets; classification was not run.',
     'service_limitation'),
    (422, 'some other validation problem', 'error'),
    ('skipped', 'no gene symbol or no c. notation', 'data'),
    ('missing', 'no result at this index', 'data'),
    ('retries_exhausted', 'no success after 4 attempts', 'transient'),
    (500, 'internal server error', 'error'),
    ('error', None, 'error'),
])
def test_categorize_error(status, message, expected):
    assert qa.categorize_error(status, message) == expected


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
# happy path and per-item errors
# ---------------------------------------------------------------------------

def test_successful_batch_returns_every_result():
    variants = [V1, V2, V4]
    client = StubAriane(lambda v: FakeResponse(200, ok_body(v)))
    errs, on_error = collect_errors()
    out = client.classify(variants, on_error)
    assert set(out) == {V1[0], V2[0], V4[0]}
    assert out[V1[0]]['predicted_label'] == 'Likely Benign'
    assert errs == []


def test_per_item_error_is_recorded_without_losing_siblings():
    """A 200 can still carry a failed item; it belongs in errors, not results."""
    variants = [V1, V2]
    body = ok_body(variants)
    body['results'][1] = {'index': 1, 'status': 'error', 'variant': 'x',
                          'error': 'could not normalise', 'result': None}
    client = StubAriane(lambda v: FakeResponse(200, body))
    errs, on_error = collect_errors()
    out = client.classify(variants, on_error)
    assert set(out) == {V1[0]}
    assert [e[0] for e in errs] == [V2[0]]
    assert 'could not normalise' in errs[0][2]


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
# 422: one bad variant must not take the batch down with it
# ---------------------------------------------------------------------------

def test_422_bisects_and_keeps_the_healthy_variants():
    """The behaviour that matters most: ARIANE rejects the whole payload when
    any one variant is unclassifiable. Bisection must isolate the offender and
    still return the other three."""
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


def test_429_backs_off_then_succeeds():
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


def test_write_summary_reads_stored_results(tmp_path):
    results = tmp_path / 'results'
    results.mkdir()
    (results / 'ga4gh_VA.one.json').write_text(json.dumps({
        '_meta': {'VRS_Digest': 'ga4gh:VA.one', 'gene': 'BRCA2',
                  'c_notation': 'c.1975T>G'},
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
    assert 'ga4gh:VA.one' in lines[1]
    assert 'Likely Benign' in lines[1]
    assert 'BS1(Strong);BP4(Supporting)' in lines[1]
    assert 'PM2' not in lines[1], 'criteria that did not apply are not "met"'
