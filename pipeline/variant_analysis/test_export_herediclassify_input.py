import json
import os

import pytest

import export_herediclassify_input as exp


# ---------------------------------------------------------------------------
# small helpers
# ---------------------------------------------------------------------------

def test_digest_to_filename_sanitizes_colon():
    assert (exp.digest_to_filename('ga4gh:VA.002K-RpU')
            == 'ga4gh_VA.002K-RpU.json')


@pytest.mark.parametrize('value,expected', [
    ('19/23', 19),
    ('19-20/23', 19),
    ('7', 7),
    (None, None),
    ('', None),
    ('abc', None),
])
def test_first_int(value, expected):
    assert exp._first_int(value) == expected


@pytest.mark.parametrize('value,expected', [
    ('0.5', True), (None, False), ('-', False),
])
def test_defined(value, expected):
    assert exp._defined(value) == expected


# ---------------------------------------------------------------------------
# build_variant_effect
# ---------------------------------------------------------------------------

def _tc(**overrides):
    tc = {
        'transcript_id': 'ENST00000357654',
        'gene_symbol': 'BRCA1',
        'hgvsc': 'ENST00000357654.9:c.5219T>G',
        'hgvsp': 'ENSP00000350283.3:p.Val1740Gly',
        'consequence_terms': ['missense_variant'],
        'exon': '19/23',
    }
    tc.update(overrides)
    return tc


def test_build_variant_effect_maps_vep_fields():
    effects = exp.build_variant_effect([_tc()], 'BRCA1')
    assert effects == [{
        'transcript': 'ENST00000357654',
        'hgvs_c': 'c.5219T>G',
        'hgvs_p': 'p.Val1740Gly',
        'variant_type': ['missense_variant'],
        'exon': 19,
        'intron': None,
    }]


def test_build_variant_effect_intron_and_null_hgvsp():
    effects = exp.build_variant_effect(
        [_tc(hgvsp=None, exon=None, intron='18/22',
             consequence_terms=['splice_donor_variant'])], 'BRCA1')
    assert effects[0]['hgvs_p'] is None
    assert effects[0]['exon'] is None
    assert effects[0]['intron'] == 18


def test_build_variant_effect_skips_non_enst_other_gene_and_missing_hgvsc():
    tcs = [
        _tc(transcript_id='NM_007294.4'),          # not Ensembl
        _tc(gene_symbol='RND1'),                   # different gene
        _tc(hgvsc=None),                           # no cDNA HGVS
        _tc(consequence_terms=[]),                 # no consequence terms
        _tc(),                                     # good
    ]
    effects = exp.build_variant_effect(tcs, 'BRCA1')
    assert len(effects) == 1


def test_build_variant_effect_empty_input():
    assert exp.build_variant_effect(None, 'BRCA1') == []


# ---------------------------------------------------------------------------
# build_gnomad
# ---------------------------------------------------------------------------

_POPULATIONS = {
    'afr': {'ac': '2', 'af': '0.0001', 'an': '41576', 'ac_hom': '1'},
    'nfe': {'ac': '10', 'af': '0.0002', 'an': '68074', 'ac_hom': '0'},
    'eas': {'ac': '0', 'af': '0.0', 'an': '5174', 'ac_hom': '2'},
}


def test_build_gnomad_absent_row_returns_none():
    assert exp.build_gnomad(None, None, None, None, None) is None


def test_build_gnomad_maps_fields_and_uppercases_subpopulation():
    block = exp.build_gnomad('0.00015', '12', '0.0001', 'nfe', _POPULATIONS)
    assert block == {
        'AF': 0.00015,
        'AC': 12,
        'AC_hom': 3,               # summed over populations
        'subpopulation': 'NFE',    # uppercased for the [A-Z]{3} schema pattern
        'popmax_AF': 0.0002,       # nfe has the highest af
        'popmax_AC': 10,
        'faf_popmax_AF': 0.0001,
    }


def test_build_gnomad_undefined_faf_and_population_fall_back():
    block = exp.build_gnomad('0.0', '0', '-', '-', _POPULATIONS)
    assert block['faf_popmax_AF'] == 0.0
    # "ALL" is HerediClassify's no-subpopulation marker and satisfies [A-Z]{3}
    assert block['subpopulation'] == 'ALL'


def test_build_gnomad_no_populations_json():
    block = exp.build_gnomad('0.0001', '5', '0.00005', 'afr', None)
    assert block['AC_hom'] == 0
    assert block['popmax_AF'] == 0.0
    assert block['popmax_AC'] == 0


# ---------------------------------------------------------------------------
# build_input_json
# ---------------------------------------------------------------------------

def _row(**overrides):
    fields = {
        'vrs_digest': 'ga4gh:VA.test1234',
        'gene': 'BRCA1',
        'hgvs_cdna': 'c.5219T>G',
        'hgvs': 'NC_000017.11:g.43057110A>C',
        'chr': '17',
        'pos': '43057110',
        'ref': 'A',
        'alt': 'C',
        'allele_frequency': '0.00015',
        'allele_count': '12',
        'faf95_popmax': '0.0001',
        'faf95_popmax_population': 'nfe',
        'populations': _POPULATIONS,
        'bayesdel': '0.31',
        'spliceai_result': '0.05',
        'prior': '0.02',
        'co_occurrence': '-',
        'segregation': None,
        'product_of_lrs': '1.5',
    }
    fields.update(overrides)
    return tuple(fields.values())


def _vep_records():
    return [{
        'most_severe_consequence': 'missense_variant',
        'transcript_consequences': [_tc()],
    }]


def test_build_input_json_full_row():
    data, reason = exp.build_input_json(_row(), _vep_records())
    assert reason is None
    assert data['_meta'] == {'VRS_Digest': 'ga4gh:VA.test1234',
                             'HGVS_cDNA': 'c.5219T>G'}
    assert data['chr'] == '17'
    assert data['pos'] == 43057110
    assert data['variant_type'] == ['missense_variant']
    assert data['pathogenicity_prediction_tools'] == {'BayesDel': 0.31}
    assert data['splicing_prediction_tools'] == {'SpliceAI': 0.05}
    assert data['gnomAD']['AF'] == 0.00015
    # exLOVD: '-' and NULL omitted, defined values kept
    assert data['prior'] == 0.02
    assert 'co-occurrence' not in data
    assert 'segregation' not in data
    assert data['multifactorial_log-likelihood'] == 1.5


def test_build_input_json_omits_gnomad_when_absent():
    row = _row(allele_frequency=None, allele_count=None, faf95_popmax=None,
               faf95_popmax_population=None, populations=None)
    data, reason = exp.build_input_json(row, _vep_records())
    assert reason is None
    assert 'gnomAD' not in data


def test_build_input_json_omits_out_of_bounds_bayesdel():
    data, _ = exp.build_input_json(_row(bayesdel='5.0'), _vep_records())
    assert 'pathogenicity_prediction_tools' not in data


def test_build_input_json_vep_error_and_empty():
    _, reason = exp.build_input_json(_row(), {'error': 'parse failed'})
    assert 'VEP error' in reason
    _, reason = exp.build_input_json(_row(), [])
    assert reason == 'no VEP records'


def test_build_input_json_no_usable_transcripts():
    records = [{'most_severe_consequence': 'missense_variant',
                'transcript_consequences': [_tc(transcript_id='NM_007294.4')]}]
    _, reason = exp.build_input_json(_row(), records)
    assert reason == 'no usable transcript consequences'


def test_build_input_json_validates_against_herediclassify_schema():
    jsonschema = pytest.importorskip('jsonschema')
    # schema_input.json is a verbatim copy of HerediClassify's
    # API/schema_input.json (GPL-3.0, commit 22905fad), kept in-repo so this
    # test needs no HerediClassify clone.
    schema_path = os.path.join(os.path.dirname(__file__), 'herediclassify',
                               'schema_input.json')
    with open(schema_path) as f:
        schema = json.load(f)

    data, _ = exp.build_input_json(_row(), _vep_records())
    jsonschema.validate(data, schema)   # raises on failure

    # gnomAD-absent variant must also validate
    row = _row(allele_frequency=None, allele_count=None, faf95_popmax=None,
               faf95_popmax_population=None, populations=None)
    data, _ = exp.build_input_json(row, _vep_records())
    jsonschema.validate(data, schema)
