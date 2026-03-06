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


def test_anw_grid_boundary_uses_matrix_gap_costs():
    """First row/column should use cumulative matrix-defined gap costs."""
    from malign.anw import nw_grids

    matrix = malign.ScoringMatrix(
        scores={
            ("A", "A"): 2.0,
            ("A", "-"): -5.0,
            ("-", "A"): -2.0,
        },
        impute_method=None,
    )

    s_grid, _ = nw_grids(["-", "A"], ["-", "A"], matrix)
    assert s_grid[0][1] == -5.0
    assert s_grid[1][0] == -2.0


# TODO: Phase 3 - Expand ANW tests
# - Test with custom matrices
# - Test edge cases (empty, single char)
# - Test scoring correctness
# - Compare with known good alignments
