"""Tests for Yen's K-Shortest Paths alignment algorithm."""

import malign


def test_yenksp_basic():
    """Test basic YenKSP alignment."""
    alms = malign.align(["ACGT", "AGCT"], k=4, method="yenksp")
    assert len(alms) <= 4
    assert len(alms[0].seqs) == 2


def test_yenksp_graph():
    """Test graph construction for YenKSP."""
    from malign.utils import DNA_MATRIX
    from malign.yenksp import compute_graph

    graph = compute_graph("ACGT", "AGCT", DNA_MATRIX)
    assert graph is not None
    assert len(graph) > 0


def test_yenksp_asymmetric_gap_orientation_regression():
    """Boundary gap edges must use the correct asymmetric symbol-gap scores."""
    from malign.yenksp import compute_graph

    matrix = malign.ScoringMatrix(
        scores={
            ("A", "A"): 2.0,
            ("A", "-"): -9.0,
            ("-", "A"): -1.0,
        },
        impute_method=None,
    )
    max_score = max(matrix.scores.values())

    # j == 0 boundary move: consume seq_a symbol against gap in seq_b.
    graph_a = compute_graph("A", "", matrix)
    assert graph_a[(0, 0)][(1, 0)] == max_score - matrix["A", "-"]

    # i == 0 boundary move: consume seq_b symbol against gap in seq_a.
    graph_b = compute_graph("", "A", matrix)
    assert graph_b[(0, 0)][(0, 1)] == max_score - matrix["-", "A"]
