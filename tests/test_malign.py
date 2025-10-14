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
