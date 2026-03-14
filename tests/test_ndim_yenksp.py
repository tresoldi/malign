"""Tests for N-dimensional YenKSP alignment."""

import malign
from malign.ndim_yenksp import build_align_nd, compute_graph_nd, ndim_yenksp_align
from malign.scoring_matrix import ScoringMatrix


def test_2seq_regression():
    """N-dim YenKSP with 2 sequences should match pairwise yenksp_align."""
    seqs = [["A", "C", "G", "T"], ["A", "G", "G", "T"]]
    matrix = ScoringMatrix.from_sequences(seqs, match=+1, gap_score=-1)

    pairwise = malign.align(seqs, k=1, method="yenksp", matrix=matrix)
    ndim = ndim_yenksp_align(seqs, k=1, matrix=matrix)

    # Both should produce valid alignments with same score
    assert len(ndim) >= 1
    assert abs(pairwise[0].score - ndim[0].score) < 0.1


def test_3seq_small():
    """3-sequence alignment produces valid results."""
    seqs = [["A", "C"], ["A", "G"], ["T", "C"]]
    matrix = ScoringMatrix.from_sequences(seqs, match=+1, gap_score=-1)

    alms = ndim_yenksp_align(seqs, k=1, matrix=matrix)
    assert len(alms) >= 1
    assert len(alms[0].seqs) == 3

    # All sequences should have same length
    lengths = {len(s) for s in alms[0].seqs}
    assert len(lengths) == 1

    # Original symbols preserved
    for i, orig in enumerate(seqs):
        aligned = [s for s in alms[0].seqs[i] if s != "-"]
        assert aligned == orig


def test_4seq_alignment():
    """4-sequence alignment."""
    seqs = [["A", "C"], ["A", "G"], ["T", "C"], ["T", "G"]]
    matrix = ScoringMatrix.from_sequences(seqs, match=+1, gap_score=-1)

    alms = ndim_yenksp_align(seqs, k=2, matrix=matrix)
    assert len(alms) >= 1
    assert len(alms[0].seqs) == 4

    for alm in alms:
        lengths = {len(s) for s in alm.seqs}
        assert len(lengths) == 1


def test_different_lengths():
    """Sequences of very different lengths."""
    seqs = [["A"], ["A", "C", "G"], ["A", "C", "G", "T", "T"]]
    matrix = ScoringMatrix.from_sequences(seqs, match=+1, gap_score=-1)

    alms = ndim_yenksp_align(seqs, k=1, matrix=matrix)
    assert len(alms) >= 1

    for i, orig in enumerate(seqs):
        aligned = [s for s in alms[0].seqs[i] if s != "-"]
        assert aligned == orig


def test_single_char_sequences():
    """Single-character sequences."""
    seqs = [["A"], ["B"], ["C"]]
    matrix = ScoringMatrix.from_sequences(seqs, match=+1, gap_score=-1)

    alms = ndim_yenksp_align(seqs, k=1, matrix=matrix)
    assert len(alms) >= 1
    assert all(len(s) >= 1 for s in alms[0].seqs)


def test_identical_sequences():
    """All-identical sequences should align perfectly."""
    seqs = [["A", "C", "G"], ["A", "C", "G"], ["A", "C", "G"]]
    matrix = ScoringMatrix.from_sequences(seqs, match=+1, gap_score=-1)

    alms = ndim_yenksp_align(seqs, k=1, matrix=matrix)
    assert len(alms) >= 1

    # Perfect diagonal alignment expected (no gaps)
    for s in alms[0].seqs:
        assert "-" not in list(s)


def test_kbest_ordering():
    """k-best alignments should have non-increasing scores."""
    seqs = [["A", "C", "G"], ["A", "G", "G"], ["T", "C", "G"]]
    matrix = ScoringMatrix.from_sequences(seqs, match=+1, gap_score=-1)

    alms = ndim_yenksp_align(seqs, k=5, matrix=matrix)
    scores = [a.score for a in alms]
    for i in range(len(scores) - 1):
        assert scores[i] >= scores[i + 1], f"Scores not sorted: {scores}"


def test_graph_node_count():
    """Graph should have expected number of nodes."""
    seqs = [["A", "C"], ["A", "G"]]
    matrix = ScoringMatrix.from_sequences(seqs, match=+1, gap_score=-1)
    padded = [[matrix.gap, *s] for s in seqs]

    graph = compute_graph_nd(padded, matrix)
    # For 2 seqs of length 2: grid is 3x3 = 9 nodes
    assert len(graph) == 9


def test_graph_3seq():
    """3-seq graph has correct structure."""
    seqs = [["A"], ["B"], ["C"]]
    matrix = ScoringMatrix.from_sequences(seqs, match=+1, gap_score=-1)
    padded = [[matrix.gap, *s] for s in seqs]

    graph = compute_graph_nd(padded, matrix)
    # Grid is 2x2x2 = 8 nodes
    assert len(graph) == 8


def test_build_align_nd():
    """build_align_nd reconstructs alignment from path."""
    seqs = [["-", "A", "C"], ["-", "B", "D"]]
    path = [(0, 0), (1, 1), (2, 2)]

    result = build_align_nd(path, seqs, "-")
    assert len(result) == 2
    assert result[0] == ["A", "C"]
    assert result[1] == ["B", "D"]


def test_build_align_nd_with_gaps():
    """build_align_nd handles gap insertions."""
    seqs = [["-", "A", "C"], ["-", "B"]]
    # Diagonal then vertical (only dim 0 advances)
    path = [(0, 0), (1, 1), (2, 1)]

    result = build_align_nd(path, seqs, "-")
    assert result[0] == ["A", "C"]
    assert result[1] == ["B", "-"]


def test_dispatch_ndim_3seq():
    """align() dispatches to ndim for 3 sequences within threshold."""
    seqs = [["A", "C", "G"], ["A", "G", "G"], ["T", "C", "G"]]
    alms = malign.align(seqs, k=1, method="anw")
    assert len(alms) >= 1
    assert len(alms[0].seqs) == 3


def test_dispatch_progressive_fallback():
    """align() falls back to progressive for sequences exceeding threshold."""
    # 5 sequences -> exceeds _NDIM_MAX_N=4
    seqs = [["A"], ["A", "C"], ["A", "C", "G"], ["A", "C", "G", "T"], ["A", "C", "G", "T", "A"]]
    alms = malign.align(seqs, k=1, method="yenksp")
    assert len(alms) >= 1
    assert len(alms[0].seqs) == 5


def test_custom_scoring_matrix():
    """N-dim with a custom asymmetric scoring matrix."""
    scores = {
        ("A", "A"): 2.0,
        ("A", "B"): -0.5,
        ("B", "A"): -1.0,  # Asymmetric
        ("B", "B"): 2.0,
        ("A", "-"): -1.0,
        ("B", "-"): -1.0,
        ("-", "A"): -1.0,
        ("-", "B"): -1.0,
    }
    matrix = ScoringMatrix(scores)
    seqs = [["A", "B"], ["B", "A"]]

    alms = ndim_yenksp_align(seqs, k=1, matrix=matrix)
    assert len(alms) >= 1
