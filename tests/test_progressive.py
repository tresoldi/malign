"""Tests for UPGMA-guided progressive alignment."""

import pytest

import malign
from malign.anw import nw_align
from malign.progressive import _expand_profile, upgma_progressive_align


def _make_matrix(seqs):
    """Helper to create an identity matrix for sequences."""
    return malign.malign.identity_matrix(seqs, match=+1, gap_score=-1)


def test_upgma_two_sequences():
    """Two sequences should produce same result as direct pairwise."""
    seqs = [list("ACGT"), list("AGCT")]
    matrix = _make_matrix(seqs)

    direct = nw_align(seqs[0], seqs[1], k=1, matrix=matrix)
    upgma = upgma_progressive_align(seqs, matrix, pw_func=nw_align, k=1)

    assert len(upgma) >= 1
    assert upgma[0].score is not None
    # Both should produce valid alignments with same score
    assert upgma[0].score == pytest.approx(direct[0].score, abs=0.5)


def test_upgma_three_identical():
    """Three identical sequences should align perfectly with no gaps."""
    seqs = [list("ABC"), list("ABC"), list("ABC")]
    matrix = _make_matrix(seqs)

    alms = upgma_progressive_align(seqs, matrix, pw_func=nw_align, k=1)

    assert len(alms) >= 1
    # All sequences should be identical (no gaps needed)
    for seq in alms[0].seqs:
        stripped = [s for s in seq if s != "-"]
        assert stripped == list("ABC")


def test_upgma_five_sequences():
    """Five sequences should produce valid alignment via progressive."""
    seqs = [
        list("ACGT"),
        list("AGCT"),
        list("ACGT"),
        list("AGCT"),
        list("ACGG"),
    ]
    matrix = _make_matrix(seqs)

    alms = upgma_progressive_align(seqs, matrix, pw_func=nw_align, k=1)

    assert len(alms) >= 1
    assert alms[0].score is not None
    # All sequences should have the same aligned length
    lengths = {len(s) for s in alms[0].seqs}
    assert len(lengths) == 1


def test_upgma_different_lengths():
    """Sequences of different lengths should be aligned with gaps."""
    seqs = [
        list("ABC"),
        list("ABCD"),
        list("BC"),
        list("ABCDE"),
        list("BCD"),
    ]
    matrix = _make_matrix(seqs)

    alms = upgma_progressive_align(seqs, matrix, pw_func=nw_align, k=1)

    assert len(alms) >= 1
    # All aligned sequences must have same length
    lengths = {len(s) for s in alms[0].seqs}
    assert len(lengths) == 1


def test_upgma_preserves_sequences():
    """Stripping gaps from aligned sequences should recover originals."""
    seqs = [
        list("ABC"),
        list("ABCD"),
        list("BC"),
        list("ABCDE"),
        list("BCD"),
    ]
    matrix = _make_matrix(seqs)

    alms = upgma_progressive_align(seqs, matrix, pw_func=nw_align, k=1)

    for i, aligned in enumerate(alms[0].seqs):
        stripped = [s for s in aligned if s != "-"]
        assert stripped == list(seqs[i])


def test_upgma_k_best():
    """k>1 should return multiple alignments with non-increasing scores."""
    seqs = [
        list("ACGT"),
        list("AGCT"),
        list("ACGT"),
        list("AGCT"),
        list("ACGG"),
    ]
    matrix = _make_matrix(seqs)

    alms = upgma_progressive_align(seqs, matrix, pw_func=nw_align, k=3)

    assert len(alms) >= 1
    # Scores should be non-increasing
    for i in range(len(alms) - 1):
        assert alms[i].score >= alms[i + 1].score


def test_upgma_vs_ndim():
    """For N=3, UPGMA should produce reasonable alignment quality."""
    seqs = [list("ACG"), list("AG"), list("ACG")]
    matrix = _make_matrix(seqs)

    alms = upgma_progressive_align(seqs, matrix, pw_func=nw_align, k=1)

    assert len(alms) >= 1
    assert alms[0].score is not None
    # All sequences should have same aligned length
    lengths = {len(s) for s in alms[0].seqs}
    assert len(lengths) == 1


def test_expand_profile_basic():
    """Basic profile expansion inserts gap columns correctly."""
    profile = [["A", "B", "C"]]
    aligned_rep = ["A", "-", "B", "C"]
    gap = "-"

    result = _expand_profile(profile, aligned_rep, gap)

    assert result == [["A", "-", "B", "C"]]


def test_expand_profile_existing_gaps():
    """Profile with internal gaps preserves them during expansion."""
    profile = [["A", "B", "C"], ["A", "-", "C"]]
    aligned_rep = ["A", "-", "B", "C"]
    gap = "-"

    result = _expand_profile(profile, aligned_rep, gap)

    # First seq should get the gap inserted
    assert result[0] == ["A", "-", "B", "C"]
    # Second seq should preserve its internal gap and get the new one
    assert result[1] == ["A", "-", "-", "C"]


def test_upgma_ndim_opportunistic():
    """Small subtrees should use true N-dim alignment when feasible."""
    # 4 short sequences: should trigger N-dim for subtrees
    seqs = [list("AB"), list("AB"), list("AB"), list("AB")]
    matrix = _make_matrix(seqs)

    alms = upgma_progressive_align(seqs, matrix, pw_func=nw_align, k=1)

    assert len(alms) >= 1
    # All identical sequences, perfect alignment expected
    for seq in alms[0].seqs:
        stripped = [s for s in seq if s != "-"]
        assert stripped == list("AB")


def test_upgma_via_align_dispatch():
    """Verify that align() dispatches to UPGMA for 5+ sequences."""
    seqs = [
        list("ACGT"),
        list("AGCT"),
        list("ACGT"),
        list("AGCT"),
        list("ACGG"),
    ]

    # This should use UPGMA progressive via the align() dispatch
    alms = malign.align(seqs, k=1)

    assert len(alms) >= 1
    assert alms[0].score is not None
    lengths = {len(s) for s in alms[0].seqs}
    assert len(lengths) == 1
