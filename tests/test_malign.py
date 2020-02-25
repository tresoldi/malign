#!/usr/bin/env python3

"""
test_malign
===========

Tests for the `malign` package.
"""

# Import Python libraries
import math
import sys
import unittest

# Impor the library itself
import malign

# TODO: use the default DNA scorer
DNA_SCORER = {
    ("A", "A"): 10,
    ("A", "G"): -1,
    ("A", "C"): -3,
    ("A", "T"): -4,
    ("G", "A"): -1,
    ("G", "G"): 7,
    ("G", "C"): -5,
    ("G", "T"): -3,
    ("C", "A"): -3,
    ("C", "G"): -5,
    ("C", "C"): 9,
    ("C", "T"): 0,
    ("T", "A"): -4,
    ("T", "G"): -3,
    ("T", "C"): 0,
    ("T", "T"): 8,
    ("A", "-"): -5,
    ("G", "-"): -5,
    ("C", "-"): -5,
    ("T", "-"): -5,
    ("-", "A"): -5,
    ("-", "G"): -5,
    ("-", "C"): -5,
    ("-", "T"): -5,
}


class TestMalign(unittest.TestCase):
    def test_dumb_align(self):
        """
        Test `dumb` pairwise alignment.
        """

        alm_a, alm_b, score = malign.align("tra", "fata", method="dumb")
        assert tuple(alm_a) == ("t", "r", "a", "-")
        assert tuple(alm_b) == ("f", "a", "t", "a")
        assert math.isclose(score, 0.75)

    def test_nw_align(self):
        """
        Test `nw` pairwise alignment.
        """

        alm_a, alm_b, score = malign.align("tra", "fata", method="nw")
        assert tuple(alm_a) == ("-", "-", "t", "r", "a")
        assert tuple(alm_b) == ("f", "a", "t", "-", "a")
        assert math.isclose(score, 0.0)

    # TODO: test class only for multiple alignment?
    def test_fill_scorer(self):
        """
        Test scorer filling.
        """

        # TODO: test from zero first
        scorer = malign.fill_scorer("AGCT", "AGCT", DNA_SCORER)
        assert math.isclose(scorer["A", "A"], DNA_SCORER["A", "A"])
        assert math.isclose(scorer["A", "T"], DNA_SCORER["A", "T"])
        assert math.isclose(scorer["A", "-"], DNA_SCORER["A", "-"])
        assert scorer["-", "-"] == 0.0

        # Remove some entries and check computation from mean values
        del scorer["A", "A"]
        del scorer["A", "T"]
        del scorer["A", "-"]

        scorer = malign.fill_scorer("AGCT", "AGCT", scorer)
        assert math.isclose(scorer["A", "A"], 6.0)
        assert math.isclose(scorer["A", "T"], -2.545454545454)
        assert math.isclose(scorer["A", "-"], -5.0)
        assert scorer["-", "-"] == 0.0

    def test_compute_graph(self):
        """
        Test scorer filling.
        """

        # TODO: test providing scorer

        dna_seq1 = [base for base in "TGGAACC"]
        dna_seq2 = [base for base in "TAGACC"]
        graph = malign.compute_graph(dna_seq1, dna_seq2, DNA_SCORER)
        assert len(graph.nodes) == 56
        assert len(graph.edges) == 122
        assert graph.edges["0:0", "1:1"]["weight"] == 2
        assert graph.edges["6:5", "7:6"]["weight"] == 1

    def test_get_aligns(self):
        """
        Test k-best alignments.
        """

        dna_seq1 = [base for base in "TGGACCCGGGAAGGTGACCCAC"]
        dna_seq2 = [base for base in "TTACCACCGGCGCGAACCCCCCCCC"]
        scorer = malign.fill_scorer("ACGT", "ACGT", DNA_SCORER)
        graph = malign.compute_graph(dna_seq1, dna_seq2, scorer)

        dest = "%i:%i" % (len(dna_seq1), len(dna_seq2))
        aligns = malign.get_aligns(graph, ("0:0", dest), dna_seq1, dna_seq2, 3)

        assert "".join(aligns[0][0][0]) == "TGG--ACC--CGGGAAGGTGACCCAC"
        assert "".join(aligns[0][0][1]) == "TTACCACCGGCGCGAACC-CCCCCCC"
        assert aligns[0][1] == 203.0

        assert "".join(aligns[1][0][0]) == "TGGAC-CCGG---GAAGGTGACCCAC"
        assert "".join(aligns[1][0][1]) == "TTACCACCGGCGCGAACC-CCCCCCC"
        assert aligns[1][1] == 204.0

        assert "".join(aligns[2][0][0]) == "TGGAC-CCGG-G--AAGGTGACCCAC"
        assert "".join(aligns[2][0][1]) == "TTACCACCGGCGCGAACC-CCCCCCC"
        assert aligns[2][1] == 205.0


if __name__ == "__main__":
    # Explicitly creating and running a test suite allows to profile it
    suite = unittest.TestLoader().loadTestsFromTestCase(TestMalign)
    unittest.TextTestRunner(verbosity=2).run(suite)
