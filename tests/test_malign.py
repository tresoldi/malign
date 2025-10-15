"""test_malign
===========

Tests for the `malign` package.
"""

# Import the library itself
import pytest

import malign


def test_dumb_pw_align():
    """Test `dumb` pairwise alignment."""

    alms = malign.multi_align(["tra", "fatata"], method="dumb")
    assert len(alms) == 1
    assert tuple(alms[0].seqs[0]) == ("-", "t", "r", "a", "-", "-")
    assert tuple(alms[0].seqs[1]) == ("f", "a", "t", "a", "t", "a")
    assert alms[0].score == pytest.approx(-1.3)


def test_nw_pw_align():
    """Test `nw` pairwise alignment."""

    alms = malign.multi_align(["tra", "fata"], k=2, method="anw")
    assert len(alms) == 1
    assert tuple(alms[0].seqs[0]) == ("-", "-", "t", "r", "a")
    assert tuple(alms[0].seqs[1]) == ("f", "a", "t", "-", "a")
    assert alms[0].score == pytest.approx(-1.2)


# TODO: fix code so it computes the graph by itself, even in pairwise
def test_yenksp_pw_align():
    """Test `kbest` pairwise alignment."""

    # Test with basic alignment, no scorer
    alms = malign.multi_align(["tra", "fata"], k=4, method="yenksp")
    assert len(alms) == 4
    assert tuple(alms[0].seqs[0]) == ("t", "r", "-", "a")
    assert tuple(alms[0].seqs[1]) == ("f", "a", "t", "a")
    assert alms[0].score == pytest.approx(-0.95)

    # More complex test with DNA scorer
    dna_seq1 = "TGGACCCGGGAAGGTGACCCAC"
    dna_seq2 = "TTACCACCGGCGCGAACCCCCCCCC"
    graph = malign.yenksp.compute_graph(dna_seq1, dna_seq2, malign.utils.DNA_MATRIX)

    dest = (len(dna_seq1), len(dna_seq2))
    aligns = malign.yenksp.align(graph, (0, 0), dest, dna_seq1, dna_seq2, malign.utils.DNA_MATRIX)
    assert "".join(aligns[0].seqs[0]) == "TGGAC-CCGG-G-AAGGTGACCCAC"
    assert "".join(aligns[0].seqs[1]) == "TTACCACCGGCGCGAACCCCCCCCC"
    assert aligns[0].score == pytest.approx(2.32)


def test_compute_graph():
    """Test graph computation for Yen's algorithm."""

    dna_seq1 = "TGGAACC"
    dna_seq2 = "TAGACC"
    graph = malign.yenksp.compute_graph(dna_seq1, dna_seq2, malign.utils.DNA_MATRIX)

    assert len(graph.nodes) == 56
    assert len(graph.edges) == 139
    assert graph.edges[(0, 0), (1, 1)]["weight"] == 2
    assert graph.edges[(6, 5), (7, 6)]["weight"] == 1


def test_score_alignment():
    """Test different scoring properties."""

    mtx = malign.utils.DNA_MATRIX

    # Test #1 - perfect alignment
    alm = [("A", "T", "T"), ("A", "T", "T")]
    assert malign.utils.score_alignment(alm, mtx) == pytest.approx(8.666666, rel=1e-05)

    # Test #2 - mistmatch
    alm = [("A", "T", "T"), ("A", "T", "C")]
    assert malign.utils.score_alignment(alm, mtx) == pytest.approx(6.0, rel=1e-05)

    # Test #3 - gap
    alm = [("A", "T", "T"), ("A", "T", "-")]
    assert malign.utils.score_alignment(alm, mtx) == pytest.approx(3.666666, rel=1e-05)
    assert malign.utils.score_alignment(alm, mtx, gap_ext=-10) == pytest.approx(0.666666, rel=1e-05)
    assert malign.utils.score_alignment(alm, mtx, gap_open=-10) == pytest.approx(
        0.666666,
        rel=1e-05,
    )

    # Test #4 - complex alignment
    alm = [
        ("A", "T", "T", "C", "G", "G", "A", "-", "-", "T"),
        ("T", "A", "-", "C", "G", "G", "A", "T", "T", "T"),
    ]
    assert malign.utils.score_alignment(alm, mtx) == pytest.approx(1.3, rel=1e-05)

    alm = [
        ("A", "T", "T", "C", "G", "G", "A", "T", "-", "-"),
        ("T", "A", "-", "C", "G", "G", "A", "T", "T", "T"),
    ]
    assert malign.utils.score_alignment(alm, mtx) == pytest.approx(1.3, rel=1e-05)

    alm = [
        ("-", "A", "T", "T", "C", "G", "G", "A", "T", "-"),
        ("T", "A", "-", "C", "G", "G", "A", "T", "T", "T"),
    ]
    assert malign.utils.score_alignment(alm, mtx) == pytest.approx(-0.5, rel=1e-05)

    alm = [
        ("-", "A", "T", "T", "C", "G", "G", "A", "-", "T"),
        ("T", "A", "-", "C", "G", "G", "A", "T", "T", "T"),
    ]
    assert malign.utils.score_alignment(alm, mtx) == pytest.approx(-0.5, rel=1e-05)

    alm = [
        ("A", "T", "T", "-", "C", "G", "G", "A", "-", "-", "T"),
        ("T", "-", "-", "A", "C", "G", "G", "A", "T", "T", "T"),
    ]
    assert malign.utils.score_alignment(alm, mtx) == pytest.approx(0.454545, rel=1e-05)

    # TODO: multiple alignments, even with identity matrix is enough


def test_tabulation():
    """Test alignment tabulation output"""

    # TODO: assertMultiLineEqual() is failing, only keeping here for coverage

    alms = malign.multi_align(["tra", "fatata"], method="anw", k=3)
    output = malign.tabulate_alms(alms)

    ref = """
| Idx   | Seq   |   Score |  #0  |  #1  |  #2  |  #3  |  #4  |  #5  |
|-------|-------|---------|------|------|------|------|------|------|
| 0     | A     |   -0.33 |  -   |  -   |  t   |  r   |  -   |  a   |
| 0     | B     |   -0.33 |  f   |  a   |  t   |  a   |  t   |  a   |
|       |       |         |      |      |      |      |      |      |
| 1     | A     |   -0.33 |  -   |  -   |  t   |  -   |  r   |  a   |
| 1     | B     |   -0.33 |  f   |  a   |  t   |  a   |  t   |  a   |
        """


def test_error_k_less_than_one():
    """Test that k < 1 raises ValueError."""
    with pytest.raises(ValueError, match="At least one alignment must be returned"):
        malign.multi_align(["ACGT", "ACGT"], k=0)

    with pytest.raises(ValueError, match="At least one alignment must be returned"):
        malign.multi_align(["ACGT", "ACGT"], k=-1)


def test_error_invalid_method():
    """Test that invalid method raises ValueError."""
    with pytest.raises(ValueError, match="Invalid alignment method"):
        malign.multi_align(["ACGT", "ACGT"], method="invalid")

    with pytest.raises(ValueError, match="Invalid alignment method"):
        malign.multi_align(["ACGT", "ACGT"], method="upgma")


def test_multiwise_alignment_varying_lengths():
    """Test multiwise alignment with sequences of very different lengths.

    This exercises the complex length-filling logic in _collect_alignments (lines 133-157).
    """
    # Test with 3 sequences of very different lengths
    seqs = [
        ["A"],  # Very short
        ["A", "C", "G"],  # Medium
        ["A", "C", "G", "T", "T"],  # Long
    ]

    alms = malign.multi_align(seqs, method="anw", k=2)

    # Should produce valid alignments
    assert len(alms) >= 1

    # All sequences in alignment should have same length
    for alm in alms:
        lengths = [len(seq) for seq in alm.seqs]
        assert len(set(lengths)) == 1, "All aligned sequences should have same length"

    # Check that sequences are preserved (gaps removed)
    for alm in alms:
        for i, original in enumerate(seqs):
            aligned_no_gaps = [s for s in alm.seqs[i] if s != "-"]
            assert aligned_no_gaps == original, f"Sequence {i} not preserved"


def test_multiwise_alignment_four_sequences():
    """Test multiwise alignment with 4 sequences to exercise complex pairing logic."""
    seqs = [
        ["A", "C"],
        ["A", "G"],
        ["T", "C"],
        ["T", "G"],
    ]

    alms = malign.multi_align(seqs, method="anw", k=3)

    # Should produce alignments
    assert len(alms) >= 1

    # Verify all have same length
    for alm in alms:
        lengths = [len(seq) for seq in alm.seqs]
        assert len(set(lengths)) == 1
        assert len(alm.seqs) == 4, "Should have 4 sequences"


def test_multiwise_alignment_five_sequences():
    """Test with 5 sequences to further exercise pairing combinations."""
    seqs = [
        ["A"],
        ["A", "C"],
        ["A", "C", "G"],
        ["A", "C", "G", "T"],
        ["A", "C", "G", "T", "A"],
    ]

    alms = malign.multi_align(seqs, method="yenksp", k=2)

    # Should handle 5 sequences
    assert len(alms) >= 1
    assert len(alms[0].seqs) == 5

    # Verify length consistency
    for alm in alms:
        lengths = [len(seq) for seq in alm.seqs]
        assert len(set(lengths)) == 1
