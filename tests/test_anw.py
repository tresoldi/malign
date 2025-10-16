"""Tests for the Asymmetric Needleman-Wunsch algorithm."""

import malign


def test_anw_basic():
    """Test basic ANW alignment."""
    alms = malign.align(["ACGT", "AGCT"], k=1, method="anw")
    assert len(alms) == 1
    assert len(alms[0].seqs) == 2


def test_anw_with_k():
    """Test ANW with multiple alignments."""
    alms = malign.align(["ACGT", "AGCT"], k=3, method="anw")
    assert len(alms) <= 3  # May be less if fewer distinct alignments exist


# TODO: Phase 3 - Expand ANW tests
# - Test with custom matrices
# - Test edge cases (empty, single char)
# - Test scoring correctness
# - Compare with known good alignments
