#!/usr/bin/env python3
# pylint: disable=no-self-use

"""
test_malign
===========

Tests for the `malign` package.
"""

# Import Python libraries
from math import isclose
import unittest

# Impor the library itself
import malign


class TestMalign(unittest.TestCase):
    """
    Suite of tests for alignment methods.
    """

    def test_dumb_pw_align(self):
        """
        Test `dumb` pairwise alignment.
        """

        alms = malign.multi_align(["tra", "fatata"], method="dumb")
        assert len(alms) == 1
        assert tuple(alms[0]["seqs"][0]) == ("-", "t", "r", "a", "-", "-")
        assert tuple(alms[0]["seqs"][1]) == ("f", "a", "t", "a", "t", "a")
        assert isclose(alms[0]["score"], 0.75)

    def test_nw_pw_align(self):
        """
        Test `nw` pairwise alignment.
        """

        alms = malign.multi_align(["tra", "fata"], k=2, method="nw")
        assert len(alms) == 1
        assert tuple(alms[0]["seqs"][0]) == ("-", "-", "t", "r", "a")
        assert tuple(alms[0]["seqs"][1]) == ("f", "a", "t", "-", "a")
        assert isclose(alms[0]["score"], -0.2)

    # TODO: fix code so it computes the graph by itself, even in pairwise
    def test_yenksp_pw_align(self):
        """
        Test `kbest` pairwise alignment.
        """

        # Test with basic alignment, no scorer
        alms = malign.multi_align(["tra", "fata"], k=4, method="yenksp")
        #        assert len(alms) == 1
        assert tuple(alms[0]["seqs"][0]) == ("-", "-", "t", "r", "a")
        assert tuple(alms[0]["seqs"][1]) == ("f", "a", "t", "-", "a")
        assert isclose(alms[0]["score"], -0.2)

        # More complex test with DNA scorer
        dna_seq1 = "TGGACCCGGGAAGGTGACCCAC"
        dna_seq2 = "TTACCACCGGCGCGAACCCCCCCCC"
        graph = malign.yenksp.compute_graph(dna_seq1, dna_seq2, malign.utils.DNA_MATRIX)

        dest = (len(dna_seq1), len(dna_seq2))
        aligns = malign.yenksp.align(graph, ((0, 0), dest), dna_seq1, dna_seq2, 3)

        assert "".join(aligns[0]["a"]) == "TGGAC-CCGG-G-AAGGTGACCCAC"
        assert "".join(aligns[0]["b"]) == "TTACCACCGGCGCGAACCCCCCCCC"
        assert isclose(aligns[0]["score_a"], 191.0)
        assert isclose(aligns[0]["score_b"], 188.0)
        assert isclose(aligns[0]["score"], 379.0)

    def test_compute_graph(self):
        """
        Test graph computation for Yen's algorithm.
        """

        dna_seq1 = "TGGAACC"
        dna_seq2 = "TAGACC"
        graph = malign.yenksp.compute_graph(dna_seq1, dna_seq2, malign.utils.DNA_MATRIX)

        assert len(graph.nodes) == 56
        assert len(graph.edges) == 139
        assert graph.edges[(0, 0), (1, 1)]["weight"] == 2
        assert graph.edges[(6, 5), (7, 6)]["weight"] == 1

    def test_score_alignment(self):
        """
        Test different scoring properties.
        """

        mtx = malign.utils.DNA_MATRIX

        # Test #1 - perfect alignment
        alm = [("A", "T", "T"), ("A", "T", "T")]
        assert isclose(malign.utils.score_alignment(alm, mtx), 26.0)

        # Test #2 - mistmatch
        alm = [("A", "T", "T"), ("A", "T", "C")]
        assert isclose(malign.utils.score_alignment(alm, mtx), 18.0)

        # Test #3 - gap
        alm = [("A", "T", "T"), ("A", "T", "-")]
        assert isclose(malign.utils.score_alignment(alm, mtx), 10.5)
        assert isclose(malign.utils.score_alignment(alm, mtx, gap_ext=-10), 1.5)
        assert isclose(malign.utils.score_alignment(alm, mtx, gap_open=-10), 2)

        # Test #4 - complex alignment
        alm = [
            ("A", "T", "T", "C", "G", "G", "A", "-", "-", "T"),
            ("T", "A", "-", "C", "G", "G", "A", "T", "T", "T"),
        ]
        assert isclose(malign.utils.score_alignment(alm, mtx), 12.0)

        alm = [
            ("A", "T", "T", "C", "G", "G", "A", "T", "-", "-"),
            ("T", "A", "-", "C", "G", "G", "A", "T", "T", "T"),
        ]
        assert isclose(malign.utils.score_alignment(alm, mtx), 12.0)

        alm = [
            ("-", "A", "T", "T", "C", "G", "G", "A", "T", "-"),
            ("T", "A", "-", "C", "G", "G", "A", "T", "T", "T"),
        ]
        assert isclose(malign.utils.score_alignment(alm, mtx), -7.5)

        alm = [
            ("-", "A", "T", "T", "C", "G", "G", "A", "-", "T"),
            ("T", "A", "-", "C", "G", "G", "A", "T", "T", "T"),
        ]
        assert isclose(malign.utils.score_alignment(alm, mtx), -7.5)

        alm = [
            ("A", "T", "T", "-", "C", "G", "G", "A", "-", "-", "T"),
            ("T", "-", "-", "A", "C", "G", "G", "A", "T", "T", "T"),
        ]
        assert isclose(malign.utils.score_alignment(alm, mtx), 2.5)

        # TODO: multiple alignments, even with identity matrix is enough


if __name__ == "__main__":
    # Explicitly creating and running a test suite allows to profile it
    suite = unittest.TestLoader().loadTestsFromTestCase(TestMalign)
    unittest.TextTestRunner(verbosity=2).run(suite)
