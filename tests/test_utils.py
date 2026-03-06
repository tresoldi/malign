"""Tests for utility functions."""

import malign


def test_tabulate_alms():
    """Test alignment tabulation."""
    alms = malign.align(["ACGT", "AGCT"], k=1)
    table = malign.tabulate_alms(alms)
    assert isinstance(table, str)
    assert len(table) > 0


def test_identity_matrix():
    """Test identity matrix creation."""
    from malign.utils import identity_matrix

    seqs = [["A", "C", "G", "T"], ["A", "C", "G", "T"]]
    matrix = identity_matrix(seqs)
    assert matrix is not None


def test_score_alignment_defaults_to_additive_objective():
    """Default scoring should match the additive matrix objective."""
    matrix = malign.ScoringMatrix(
        scores={
            ("A", "A"): 2.0,
            ("A", "-"): -1.0,
            ("-", "A"): -1.0,
        },
        impute_method=None,
    )

    from malign.utils import score_alignment

    score = score_alignment([["A"], ["A"]], matrix)
    assert score == 2.0


def test_score_alignment_optional_normalize_and_affine_penalties():
    """Optional flags keep normalized/penalized scoring available."""
    matrix = malign.ScoringMatrix(
        scores={
            ("A", "A"): 2.0,
            ("A", "-"): -1.0,
            ("-", "A"): -1.0,
        },
        impute_method=None,
    )

    from malign.utils import score_alignment

    score = score_alignment(
        [["A", "-"], ["A", "A"]],
        matrix,
        normalize=True,
        gap_ext=-1.0,
        gap_open=-1.0,
    )
    assert score == -0.5
