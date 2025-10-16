"""Tests for the dumb baseline alignment method."""

import malign


def test_dumb_basic():
    """Test basic dumb alignment."""
    alms = malign.align(["ACGT", "AGCT"], k=1, method="dumb")
    assert len(alms) == 1
    assert len(alms[0].seqs) == 2


# TODO: Phase 3 - Add more dumb method tests
# - Verify it's always worse than ANW
# - Use as baseline for benchmarking
