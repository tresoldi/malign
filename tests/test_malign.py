#!/usr/bin/env python3
# pylint: disable=no-self-use

"""
test_malign
===========

Tests for the `malign` package.
"""

# Import Python libraries
import math
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
        assert math.isclose(alms[0]["score"], 0.75)

    def test_nw_pw_align(self):
        """
        Test `nw` pairwise alignment.
        """

        alms = malign.multi_align(["tra", "fata"], k=2, method="nw")
        assert len(alms) == 1
        assert tuple(alms[0]["seqs"][0]) == ("-", "-", "t", "r", "a")
        assert tuple(alms[0]["seqs"][1]) == ("f", "a", "t", "-", "a")
        assert math.isclose(alms[0]["score"], -0.2)

    # TODO: rename to yenksp
    # TODO: fix code so it computes the graph by itself, even in pairwise
    def test_kbest_pw_align(self):
        """
        Test `kbest` pairwise alignment.
        """

        # Test with basic alignment, no scorer
        alms = malign.multi_align(["tra", "fata"], k=4, method="yenksp")
        #        assert len(alms) == 1
        assert tuple(alms[0]["seqs"][0]) == ("-", "-", "t", "r", "a")
        assert tuple(alms[0]["seqs"][1]) == ("f", "a", "t", "-", "a")
        assert math.isclose(alms[0]["score"], -0.2)

        # More complex test with DNA scorer
        dna_seq1 = "TGGACCCGGGAAGGTGACCCAC"
        dna_seq2 = "TTACCACCGGCGCGAACCCCCCCCC"
        graph = malign.kbest.compute_graph(dna_seq1, dna_seq2, malign.utils.DNA_MATRIX)

        dest = "%i:%i" % (len(dna_seq1), len(dna_seq2))
        aligns = malign.kbest.align(graph, ("0:0", dest), dna_seq1, dna_seq2, 3)

        assert "".join(aligns[0]["a"]) == "TGGAC-CCGG-G-AAGGTGACCCAC"
        assert "".join(aligns[0]["b"]) == "TTACCACCGGCGCGAACCCCCCCCC"
        assert math.isclose(aligns[0]["score_a"], 191.0)
        assert math.isclose(aligns[0]["score_b"], 188.0)
        assert math.isclose(aligns[0]["score"], 379.0)

    def test_compute_graph(self):
        """
        Test graph computation for Yen's algorithm.
        """

        dna_seq1 = "TGGAACC"
        dna_seq2 = "TAGACC"
        graph = malign.kbest.compute_graph(dna_seq1, dna_seq2, malign.utils.DNA_MATRIX)

        assert len(graph.nodes) == 56
        assert len(graph.edges) == 139
        assert graph.edges["0:0", "1:1"]["weight"] == 2
        assert graph.edges["6:5", "7:6"]["weight"] == 1


if __name__ == "__main__":
    # Explicitly creating and running a test suite allows to profile it
    suite = unittest.TestLoader().loadTestsFromTestCase(TestMalign)
    unittest.TextTestRunner(verbosity=2).run(suite)
