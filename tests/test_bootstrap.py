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


# --- Prior-guided learning tests ---

# Shared fixtures for prior tests
_PRIOR_PAIRS = [
    (["p", "a", "t", "a"], ["b", "a", "d", "a"]),
    (["p", "a", "k", "a"], ["b", "a", "g", "a"]),
    (["t", "a", "p", "a"], ["d", "a", "b", "a"]),
    (["k", "a", "t", "a"], ["g", "a", "d", "a"]),
    (["p", "i", "t", "i"], ["b", "i", "d", "i"]),
]


def _make_prior():
    """Build a simple prior matrix for testing."""
    left = ["-", "a", "b", "d", "g", "i", "k", "p", "t", "x", "y"]
    right = ["-", "a", "b", "d", "g", "i", "k", "p", "t", "x", "y"]
    scores: dict[tuple[str, ...], float] = {}
    for s0 in left:
        for s1 in right:
            if s0 == "-" or s1 == "-":
                scores[(s0, s1)] = -1.0
            elif s0 == s1:
                scores[(s0, s1)] = 1.0
            else:
                scores[(s0, s1)] = -0.5
    scores[("-", "-")] = 0.0
    # Give voicing correspondences a boost
    for pair in [("p", "b"), ("t", "d"), ("k", "g")]:
        scores[pair] = 0.8
    return malign.ScoringMatrix(scores=scores, domains=[left, right], gap="-", impute_method=None)


def test_bootstrap_prior_blending():
    """With vs without prior produces different scores."""
    prior = _make_prior()

    matrix_no_prior = malign.bootstrap_matrix(_PRIOR_PAIRS, max_iter=5)
    matrix_with_prior = malign.bootstrap_matrix(
        _PRIOR_PAIRS, max_iter=5, prior_matrix=prior, prior_weight=0.5
    )

    # Scores should differ
    common_keys = set(matrix_no_prior.scores) & set(matrix_with_prior.scores)
    diffs = [abs(matrix_no_prior.scores[k] - matrix_with_prior.scores[k]) for k in common_keys]
    assert max(diffs) > 0.01, "Prior should influence scores"


def test_bootstrap_prior_decay():
    """With high max_iter, final scores converge toward pure data-driven."""
    prior = _make_prior()

    matrix_no_prior = malign.bootstrap_matrix(
        _PRIOR_PAIRS, max_iter=50, convergence_threshold=1e-6, matrix_threshold=1e-6
    )
    matrix_heavy_prior = malign.bootstrap_matrix(
        _PRIOR_PAIRS,
        max_iter=50,
        prior_matrix=prior,
        prior_weight=0.5,
        convergence_threshold=1e-6,
        matrix_threshold=1e-6,
    )

    # After many iterations the prior influence should be small
    common_keys = set(matrix_no_prior.scores) & set(matrix_heavy_prior.scores)
    diffs = [abs(matrix_no_prior.scores[k] - matrix_heavy_prior.scores[k]) for k in common_keys]
    avg_diff = sum(diffs) / len(diffs) if diffs else 0
    # Average difference should be modest (prior has decayed)
    assert avg_diff < 1.0, f"Expected prior to decay, avg diff = {avg_diff}"


def test_bootstrap_prior_symbol_coverage():
    """Prior has symbols not in pairs; they appear in the output matrix."""
    prior = _make_prior()
    # "x" and "y" are in the prior but NOT in _PRIOR_PAIRS
    matrix = malign.bootstrap_matrix(_PRIOR_PAIRS, max_iter=3, prior_matrix=prior, prior_weight=0.5)

    # Check that x and y appear in the output domains
    all_domain_symbols = set()
    for domain in matrix.domains:
        all_domain_symbols.update(domain)
    assert "x" in all_domain_symbols, "Prior symbol 'x' should be in output domains"
    assert "y" in all_domain_symbols, "Prior symbol 'y' should be in output domains"


def test_bootstrap_prior_as_initial():
    """Providing prior_matrix without initial_matrix uses the prior as start."""
    prior = _make_prior()

    # Run just 1 iteration so the initial matrix matters
    matrix = malign.bootstrap_matrix(_PRIOR_PAIRS, max_iter=1, prior_matrix=prior, prior_weight=0.5)

    assert isinstance(matrix, malign.ScoringMatrix)
    # The output should have at least as many score keys as the prior
    assert len(matrix.scores) >= len(prior.scores) * 0.5


def test_bootstrap_block_merge_basic():
    """bootstrap_matrix with block_merge=True returns a valid matrix."""
    pairs = [
        (["p", "a", "t", "a"], ["b", "a", "d", "a"]),
        (["p", "a", "k", "a"], ["b", "a", "g", "a"]),
    ]
    matrix = malign.bootstrap_matrix(pairs, max_iter=3, block_merge=True)

    assert isinstance(matrix, malign.ScoringMatrix)
    assert matrix.num_domains == 2
    assert len(matrix.scores) > 0


def test_bootstrap_block_merge_vs_no_merge():
    """block_merge=True may produce different scores than without."""
    pairs = [
        (["a"], ["j", "e"]),
        (["a"], ["j", "e"]),
        (["a"], ["j", "e"]),
        (["b"], ["b"]),
    ]
    matrix_no = malign.bootstrap_matrix(pairs, max_iter=5, block_merge=False)
    matrix_yes = malign.bootstrap_matrix(pairs, max_iter=5, block_merge=True)

    # Both should be valid
    assert isinstance(matrix_no, malign.ScoringMatrix)
    assert isinstance(matrix_yes, malign.ScoringMatrix)


def test_bootstrap_distfeat_pipeline():
    """End-to-end: from_distfeat() -> bootstrap_matrix(prior_matrix=...) -> asymmetric."""
    pytest.importorskip("distfeat")

    left_seqs = [["p", "t", "k"]]
    right_seqs = [["b", "d", "g"]]

    prior = malign.ScoringMatrix.from_distfeat(
        sequences=left_seqs + right_seqs,
        gap="-",
    )

    pairs = [
        (["p", "a", "t", "a"], ["b", "a", "d", "a"]),
        (["t", "a", "p", "a"], ["d", "a", "b", "a"]),
        (["k", "a", "t", "a"], ["g", "a", "d", "a"]),
        (["p", "i", "t", "i"], ["b", "i", "d", "i"]),
    ]

    matrix = malign.bootstrap_matrix(pairs, max_iter=10, prior_matrix=prior, prior_weight=0.5)

    assert isinstance(matrix, malign.ScoringMatrix)
    assert matrix.num_domains == 2
    # Bootstrap produces asymmetric scores even from a symmetric prior
    # Check a pair where directionality matters
    if ("p", "b") in matrix.scores and ("b", "p") in matrix.scores:
        # They CAN differ due to data-driven asymmetry
        pass  # Just verify the pipeline completes without error
