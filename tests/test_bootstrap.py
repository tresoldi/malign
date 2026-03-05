"""Tests for unsupervised bootstrap matrix learning."""

import pytest

import malign


def test_bootstrap_basic():
    """bootstrap_matrix returns a valid ScoringMatrix."""
    pairs = [
        ("ABC", "ABC"),
        ("ABCD", "ABCE"),
    ]
    matrix = malign.bootstrap_matrix(pairs, max_iter=3)

    assert isinstance(matrix, malign.ScoringMatrix)
    assert matrix.num_domains == 2
    assert matrix.gap == "-"
    assert len(matrix.scores) > 0


def test_bootstrap_convergence():
    """bootstrap_matrix converges before max_iter with sufficient data."""
    pairs = [
        ("ABC", "ABC"),
        ("ABC", "ABC"),
        ("DEF", "DEF"),
        ("DEF", "DEF"),
    ]
    matrix = malign.bootstrap_matrix(pairs, max_iter=50, convergence_threshold=0.01)

    assert isinstance(matrix, malign.ScoringMatrix)


def test_bootstrap_preserves_asymmetry():
    """Frequent (A,B) vs rare (B,A) should produce asymmetric scores."""
    pairs = [
        ("AB", "BB"),
        ("AB", "BB"),
        ("AB", "BB"),
        ("AB", "BB"),
        ("AB", "BB"),
        ("BA", "AA"),
    ]
    matrix = malign.bootstrap_matrix(pairs, max_iter=10)

    # Scores for (A,B) and (B,A) should differ
    assert ("A", "B") in matrix.scores
    assert ("B", "A") in matrix.scores
    assert matrix.scores[("A", "B")] != matrix.scores[("B", "A")]


def test_bootstrap_learned_matrix_aligns():
    """Learned matrix should produce valid alignments."""
    pairs = [
        ("ABC", "ABC"),
        ("ABCD", "ABCE"),
        ("BC", "BD"),
    ]
    matrix = malign.bootstrap_matrix(pairs, max_iter=5)

    # Use the learned matrix for alignment
    alms = malign.align(["ABC", "ABD"], k=1, matrix=matrix)

    assert len(alms) >= 1
    assert alms[0].score is not None


def test_bootstrap_empty_pairs():
    """Empty pairs list should raise ValueError."""
    with pytest.raises(ValueError, match="At least one pair"):
        malign.bootstrap_matrix([])


def test_bootstrap_single_pair():
    """Single pair should work without errors."""
    pairs = [("AB", "CD")]
    matrix = malign.bootstrap_matrix(pairs, max_iter=3)

    assert isinstance(matrix, malign.ScoringMatrix)
    assert matrix.num_domains == 2


def test_bootstrap_custom_gap():
    """Custom gap symbol should be respected."""
    pairs = [
        ("ABC", "ABC"),
        ("AB", "AC"),
    ]
    matrix = malign.bootstrap_matrix(pairs, max_iter=3, gap=".")

    assert matrix.gap == "."
    assert (".", ".") in matrix.scores
    assert matrix.scores[(".", ".")] == 0.0


def test_bootstrap_integration():
    """Full pipeline: generate pairs, learn matrix, align."""
    # Simulate phonological sound changes: p->b, t->d (voicing)
    pairs = [
        (["p", "a", "t", "a"], ["b", "a", "d", "a"]),
        (["p", "a", "k", "a"], ["b", "a", "g", "a"]),
        (["t", "a", "p", "a"], ["d", "a", "b", "a"]),
        (["k", "a", "t", "a"], ["g", "a", "d", "a"]),
        (["p", "i", "t", "i"], ["b", "i", "d", "i"]),
    ]

    matrix = malign.bootstrap_matrix(pairs, max_iter=10)

    # Use learned matrix to align a new pair
    alms = malign.align(
        [["p", "a", "t"], ["b", "a", "d"]],
        k=1,
        matrix=matrix,
    )

    assert len(alms) >= 1
    assert alms[0].score is not None

    # p->b should have a positive log-odds score (frequently co-occur)
    assert matrix.scores.get(("p", "b"), 0) > matrix.scores.get(("p", "g"), float("-inf"))
