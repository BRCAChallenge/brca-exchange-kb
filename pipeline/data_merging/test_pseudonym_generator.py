"""
Tests for pseudonym_generator HGVS_cDNA fill-in behaviour.

Focus: ensure_mane_transcript_cdna and the surrounding process_rows flow,
specifically the cases where HGVS_cDNA is not supplied by the source VCF
and must be derived from SeqRepo / ClinGen.
"""

import copy
from unittest.mock import MagicMock, patch, call
import pytest
import hgvs.exceptions

from data_merging.pseudonym_generator import (
    ensure_mane_transcript_cdna,
    HGVS_CDNA_COL,
    PYHGVS_CDNA_COL,
    PYHGVS_GENOMIC_COORDINATE_38_COL,
    REFERENCE_SEQUENCE_COL,
)

MANE = 'NM_007294.4'


def _row(hgvs_cdna='-', pyhgvs_cdna=None, ref_seq='-',
         genomic_38='NC_000017.11:g.43094692G>A'):
    """Return a minimal row dict as produced by _load_rows_from_db / process_rows."""
    return {
        HGVS_CDNA_COL:                   hgvs_cdna,
        PYHGVS_CDNA_COL:                  pyhgvs_cdna,
        REFERENCE_SEQUENCE_COL:           ref_seq,
        PYHGVS_GENOMIC_COORDINATE_38_COL: genomic_38,
    }


def _mock_hgvs_proc(parsed_obj=None):
    proc = MagicMock()
    proc.hgvs_parser.parse.return_value = parsed_obj or MagicMock()
    return proc


# ---------------------------------------------------------------------------
# ensure_mane_transcript_cdna
# ---------------------------------------------------------------------------

class TestEnsureManeTranscriptCdna:

    def test_seqrepo_fills_missing_hgvs(self):
        """When HGVS_cDNA is '-' and SeqRepo succeeds, the row gets the cDNA value."""
        row = _row(hgvs_cdna='-', pyhgvs_cdna=None)
        hgvs_proc = _mock_hgvs_proc()
        normalizer = MagicMock()
        am38 = MagicMock()

        cdna_result = MagicMock()
        cdna_result.__str__ = lambda self: f'{MANE}:c.5266dupC'
        am38.g_to_c.return_value = cdna_result
        normalizer.normalize.return_value = cdna_result

        ensure_mane_transcript_cdna(row, MANE, hgvs_proc, normalizer, am38)

        assert row[HGVS_CDNA_COL] == 'c.5266dupC'
        assert row[REFERENCE_SEQUENCE_COL] == MANE

    def test_seqrepo_fills_missing_hgvs_brca2(self):
        """Same check for a BRCA2 MANE transcript."""
        mane2 = 'NM_000059.4'
        row = _row(hgvs_cdna='-', pyhgvs_cdna=None,
                   genomic_38='NC_000013.11:g.32340300A>G')
        hgvs_proc = _mock_hgvs_proc()
        normalizer = MagicMock()
        am38 = MagicMock()

        cdna_result = MagicMock()
        cdna_result.__str__ = lambda self: f'{mane2}:c.1114A>G'
        am38.g_to_c.return_value = cdna_result
        normalizer.normalize.return_value = cdna_result

        ensure_mane_transcript_cdna(row, mane2, hgvs_proc, normalizer, am38)

        assert row[HGVS_CDNA_COL] == 'c.1114A>G'

    def test_clingen_cdna_accepted_without_seqrepo_call(self):
        """When ClinGen already returned the MANE cDNA, SeqRepo is not called."""
        pyhgvs_cdna = f'{MANE}:c.5266dupC'
        row = _row(hgvs_cdna='-', pyhgvs_cdna=pyhgvs_cdna, ref_seq=MANE)
        hgvs_proc = _mock_hgvs_proc()
        normalizer = MagicMock()
        am38 = MagicMock()

        ensure_mane_transcript_cdna(row, MANE, hgvs_proc, normalizer, am38)

        # SeqRepo must not have been called; value comes from ClinGen result
        am38.g_to_c.assert_not_called()
        assert row[HGVS_CDNA_COL] == 'c.5266dupC'

    def test_seqrepo_failure_preserves_existing_source_hgvs(self):
        """When SeqRepo fails but the row already has a VCF-sourced HGVS_cDNA,
        the existing value is kept (the key bug fix)."""
        existing = 'c.-40+1G>A'
        row = _row(hgvs_cdna=existing, pyhgvs_cdna=None)
        hgvs_proc = _mock_hgvs_proc()
        normalizer = MagicMock()
        am38 = MagicMock()
        am38.g_to_c.side_effect = hgvs.exceptions.HGVSInvalidIntervalError('out of range')

        ensure_mane_transcript_cdna(row, MANE, hgvs_proc, normalizer, am38)

        assert row[HGVS_CDNA_COL] == existing, (
            "SeqRepo failure must not overwrite a pre-existing HGVS_cDNA value"
        )

    def test_seqrepo_failure_no_existing_hgvs_stays_dash(self):
        """When SeqRepo fails and there was no pre-existing HGVS_cDNA,
        the value remains '-' (correct for gnomAD-only deep intronic variants)."""
        row = _row(hgvs_cdna='-', pyhgvs_cdna=None)
        hgvs_proc = _mock_hgvs_proc()
        normalizer = MagicMock()
        am38 = MagicMock()
        am38.g_to_c.side_effect = hgvs.exceptions.HGVSInvalidIntervalError('out of range')

        ensure_mane_transcript_cdna(row, MANE, hgvs_proc, normalizer, am38)

        assert row[HGVS_CDNA_COL] == '-'

    def test_seqrepo_unsupported_operation_preserves_existing(self):
        """HGVSUnsupportedOperationError is also caught and handled correctly."""
        existing = 'c.67+58A>C'
        row = _row(hgvs_cdna=existing, pyhgvs_cdna=None)
        hgvs_proc = _mock_hgvs_proc()
        normalizer = MagicMock()
        am38 = MagicMock()
        am38.g_to_c.side_effect = hgvs.exceptions.HGVSUnsupportedOperationError('unsupported')

        ensure_mane_transcript_cdna(row, MANE, hgvs_proc, normalizer, am38)

        assert row[HGVS_CDNA_COL] == existing

    def test_seqrepo_data_not_available_preserves_existing(self):
        """HGVSDataNotAvailableError is also caught and handled correctly."""
        existing = 'c.67+82C>G'
        row = _row(hgvs_cdna=existing, pyhgvs_cdna=None)
        hgvs_proc = _mock_hgvs_proc()
        normalizer = MagicMock()
        am38 = MagicMock()
        am38.g_to_c.side_effect = hgvs.exceptions.HGVSDataNotAvailableError('no data')

        ensure_mane_transcript_cdna(row, MANE, hgvs_proc, normalizer, am38)

        assert row[HGVS_CDNA_COL] == existing

    def test_seqrepo_overrides_non_mane_pyhgvs_cdna(self):
        """When ClinGen returns a cDNA on a non-MANE transcript, SeqRepo remaps it."""
        non_mane = 'NM_007294.3:c.5266dupC'
        row = _row(hgvs_cdna='-', pyhgvs_cdna=non_mane)
        hgvs_proc = _mock_hgvs_proc()
        normalizer = MagicMock()
        am38 = MagicMock()

        cdna_result = MagicMock()
        cdna_result.__str__ = lambda self: f'{MANE}:c.5266dupC'
        am38.g_to_c.return_value = cdna_result
        normalizer.normalize.return_value = cdna_result

        ensure_mane_transcript_cdna(row, MANE, hgvs_proc, normalizer, am38)

        am38.g_to_c.assert_called_once()
        assert row[HGVS_CDNA_COL] == 'c.5266dupC'
        assert row[REFERENCE_SEQUENCE_COL] == MANE
