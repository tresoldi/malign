"""Tests for the Alignment class and validation metrics."""


def test_alignment_creation():
    """Test basic Alignment object creation."""
    from malign.alignment import Alignment

    alm = Alignment(seqs=[["A", "T"], ["A", "C"]], score=1.0)
    assert len(alm) == 2
    assert alm.score == 1.0


def test_alignment_indexing():
    """Test Alignment indexing."""
    from malign.alignment import Alignment

    alm = Alignment(seqs=[["A", "T"], ["A", "C"]], score=1.0)
    assert alm[0] == ["A", "T"]
    assert alm[1] == ["A", "C"]


# TODO: Phase 2 - Add tests for validation metrics
# def test_sum_of_pairs_score():
#     pass
#
# def test_entropy():
#     pass
#
# def test_consensus():
#     pass
#
# def test_gap_statistics():
#     pass
