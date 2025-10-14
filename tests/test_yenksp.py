"""Tests for Yen's K-Shortest Paths alignment algorithm."""

import malign


def test_yenksp_basic():
    """Test basic YenKSP alignment."""
    alms = malign.multi_align(["ACGT", "AGCT"], k=4, method="yenksp")
    assert len(alms) <= 4
    assert len(alms[0].seqs) == 2


def test_yenksp_graph():
    """Test graph construction for YenKSP."""
    from malign.utils import DNA_MATRIX
    from malign.yenksp import compute_graph

    graph = compute_graph("ACGT", "AGCT", DNA_MATRIX)
    assert graph is not None
    assert len(graph.nodes) > 0


# TODO: Phase 3 - Expand YenKSP tests
# - Test with various k values
# - Compare diversity vs ANW
# - Test graph properties
# - Performance benchmarks
