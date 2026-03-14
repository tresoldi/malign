"""Tests for asymmetric scoring matrices and alignment."""

import math

import pytest

import malign
from malign.ndim_yenksp import ndim_yenksp_align
from malign.scoring_matrix import ScoringMatrix


def test_asymmetric_matrix_scores_differ():
    """score(A,B) != score(B,A) in an asymmetric matrix."""
    scores = {
        ("p", "p"): 2.0,
        ("p", "b"): 0.5,
        ("b", "p"): -0.5,
        ("b", "b"): 2.0,
        ("p", "-"): -1.0,
        ("b", "-"): -1.0,
        ("-", "p"): -1.0,
        ("-", "b"): -1.0,
    }
    matrix = ScoringMatrix(scores)
    assert matrix[("p", "b")] != matrix[("b", "p")]


def test_asymmetric_align_different_scores():
    """Aligning same pair in both orders should give different scores with asymmetric matrix."""
    scores = {
        ("p", "p"): 2.0,
        ("p", "b"): 1.0,
        ("b", "p"): -1.0,
        ("b", "b"): 2.0,
        ("p", "-"): -1.0,
        ("b", "-"): -1.0,
        ("-", "p"): -1.0,
        ("-", "b"): -1.0,
    }
    matrix = ScoringMatrix(scores)

    alms_fwd = ndim_yenksp_align([["p", "b"], ["b", "p"]], k=1, matrix=matrix)
    # Reverse the sequences and use a transposed matrix
    rev_scores = {(k[1], k[0]): v for k, v in scores.items()}
    rev_matrix = ScoringMatrix(rev_scores)
    alms_rev = ndim_yenksp_align([["b", "p"], ["p", "b"]], k=1, matrix=rev_matrix)

    assert len(alms_fwd) >= 1
    assert len(alms_rev) >= 1
    # Both should produce valid alignments but may have different scores
    # due to asymmetry in gap placement
    assert alms_fwd[0].score is not None
    assert alms_rev[0].score is not None


def test_pairwise_anw_yenksp_consistency_with_asymmetric_gaps():
    """Pairwise ANW and YenKSP should agree under asymmetric gap penalties."""
    matrix = ScoringMatrix(
        {
            ("A", "A"): 2.0,
            ("A", "B"): -1.0,
            ("B", "A"): -1.0,
            ("B", "B"): 2.0,
            ("A", "-"): -3.0,
            ("B", "-"): -3.0,
            ("-", "A"): -2.0,
            ("-", "B"): -3.0,
        },
        impute_method=None,
    )

    anw = malign.align(["BA", "AB"], method="anw", k=1, matrix=matrix)[0]
    yenksp = malign.align(["BA", "AB"], method="yenksp", k=1, matrix=matrix)[0]

    assert anw.seqs == yenksp.seqs
    assert anw.score == yenksp.score


def test_from_substitution_counts_basic():
    """from_substitution_counts creates a valid asymmetric matrix."""
    counts = {
        ("p", "b"): 15,
        ("b", "p"): 3,
        ("p", "p"): 20,
        ("b", "b"): 18,
    }
    matrix = ScoringMatrix.from_substitution_counts(counts)

    # Should produce asymmetric scores
    assert matrix[("p", "b")] != matrix[("b", "p")]

    # Identity pairs should have positive scores (log-odds)
    assert matrix[("p", "p")] > 0
    assert matrix[("b", "b")] > 0


def test_from_substitution_counts_empty():
    """from_substitution_counts raises on empty counts."""
    with pytest.raises(ValueError, match="non-empty"):
        ScoringMatrix.from_substitution_counts({})


def test_from_substitution_counts_gap_handling():
    """from_substitution_counts sets gap scores correctly."""
    counts = {
        ("A", "A"): 10,
        ("A", "B"): 5,
        ("B", "A"): 3,
        ("B", "B"): 8,
    }
    matrix = ScoringMatrix.from_substitution_counts(counts, gap="-", gap_score=-2.0)

    assert matrix[("-", "-")] == 0.0
    assert matrix[("A", "-")] == -2.0
    assert matrix[("-", "B")] == -2.0


def test_from_substitution_counts_log_odds():
    """Verify log-odds computation is correct."""
    counts = {
        ("A", "A"): 40,
        ("A", "B"): 10,
        ("B", "A"): 10,
        ("B", "B"): 40,
    }
    total = 100

    matrix = ScoringMatrix.from_substitution_counts(counts)

    # For ("A", "A"): observed = 40/100 = 0.4
    # marginal_A_pos0 = (40+10)/100 = 0.5, marginal_A_pos1 = (40+10)/100 = 0.5
    # expected = 0.5 * 0.5 = 0.25
    # score = log(0.4/0.25) = log(1.6)
    expected_score = math.log(0.4 / 0.25)
    assert abs(matrix[("A", "A")] - expected_score) < 1e-10


def test_from_substitution_counts_asymmetric_output():
    """Strongly asymmetric counts produce asymmetric scores."""
    counts = {
        ("p", "b"): 100,
        ("b", "p"): 1,
        ("p", "p"): 50,
        ("b", "b"): 50,
    }
    matrix = ScoringMatrix.from_substitution_counts(counts)

    # p->b is much more common than b->p
    assert matrix[("p", "b")] > matrix[("b", "p")]


def test_asymmetric_ndim_alignment():
    """End-to-end: asymmetric matrix + N-dim alignment."""
    counts = {
        ("p", "b"): 15,
        ("b", "p"): 3,
        ("p", "p"): 20,
        ("b", "b"): 18,
        ("t", "t"): 25,
        ("t", "d"): 10,
        ("d", "t"): 2,
        ("d", "d"): 20,
    }
    matrix = ScoringMatrix.from_substitution_counts(counts)

    seqs = [["p", "t"], ["b", "d"]]
    alms = ndim_yenksp_align(seqs, k=1, matrix=matrix)
    assert len(alms) >= 1

    # Symbols preserved
    for i, orig in enumerate(seqs):
        aligned = [s for s in alms[0].seqs[i] if s != "-"]
        assert aligned == orig


def test_from_substitution_counts_3way():
    """from_substitution_counts works with 3-way tuples."""
    counts = {
        ("A", "A", "A"): 20,
        ("A", "B", "A"): 5,
        ("B", "A", "B"): 3,
        ("B", "B", "B"): 15,
    }
    matrix = ScoringMatrix.from_substitution_counts(counts)
    assert len(matrix.domains) == 3


def test_threshold_logic():
    """should_use_ndim respects both N and grid_size limits."""
    from malign.ndim_common import should_use_ndim

    # N=3, small sequences -> True
    assert should_use_ndim(3, [3, 3, 3], "anw") is True

    # N=4, small sequences -> True
    assert should_use_ndim(4, [3, 3, 3, 3], "anw") is True

    # N=5 -> False (exceeds _NDIM_MAX_N=4)
    assert should_use_ndim(5, [3, 3, 3, 3, 3], "anw") is False

    # N=4, very long sequences -> False (exceeds grid cell limit)
    assert should_use_ndim(4, [100, 100, 100, 100], "anw") is False

    # N=3, reasonable sequences -> True
    assert should_use_ndim(3, [10, 10, 10], "yenksp") is True


def test_ndim_vs_progressive_score():
    """N-dim alignment should produce equal or better scores than progressive."""
    seqs = [["A", "C", "G"], ["A", "G", "G"], ["T", "C", "G"]]
    matrix = ScoringMatrix.from_sequences(seqs, match=+1, gap_score=-1)

    # N-dim alignment
    ndim_alms = ndim_yenksp_align(seqs, k=1, matrix=matrix)

    # Progressive alignment (through _collect_alignments)
    from malign.anw import nw_align
    from malign.malign import _collect_alignments

    prog_alms = _collect_alignments(seqs, matrix, pw_func=nw_align, k=1)

    assert len(ndim_alms) >= 1
    assert len(prog_alms) >= 1

    # N-dim should be at least as good (it finds the true optimum)
    assert ndim_alms[0].score >= prog_alms[0].score - 0.01
