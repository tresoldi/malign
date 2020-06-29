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


class TestMalign(unittest.TestCase):
    def test_dumb_pw_align(self):
        """
        Test `dumb` pairwise alignment.
        """

        alms = malign.pw_align("tra", "fatata", method="dumb")
        assert len(alms) == 1
        assert tuple(alms[0]["a"]) == ("-", "t", "r", "a", "-", "-")
        assert tuple(alms[0]["b"]) == ("f", "a", "t", "a", "t", "a")
        assert math.isclose(alms[0]["score"], 0.5)

    def test_nw_pw_align(self):
        """
        Test `nw` pairwise alignment.
        """

        alms = malign.pw_align("tra", "fata", k=2, method="nw")
        assert len(alms) == 1
        assert tuple(alms[0]["a"]) == ("-", "-", "t", "r", "a")
        assert tuple(alms[0]["b"]) == ("f", "a", "t", "-", "a")
        assert math.isclose(alms[0]["score"], -3.5)

    def test_kbest_pw_align(self):
        """
        Test `kbest` pairwise alignment.
        """

        # Test with basic alignment, no scorer
        alms = malign.pw_align("tra", "fata", k=4, method="kbest")
        assert len(alms) == 4
        assert tuple(alms[0]["a"]) == ("-", "t", "r", "a")
        assert tuple(alms[0]["b"]) == ("f", "a", "t", "a")
        assert math.isclose(alms[0]["score"], 5.0)

        # More complex test with DNA scorer
        dna_seq1 = [base for base in "TGGACCCGGGAAGGTGACCCAC"]
        dna_seq2 = [base for base in "TTACCACCGGCGCGAACCCCCCCCC"]
        scorer = malign.utils.fill_matrix("ACGT", "ACGT", malign.utils.DNA_MATRIX)
        graph = malign.kbest.compute_graph(dna_seq1, dna_seq2, scorer)

        dest = "%i:%i" % (len(dna_seq1), len(dna_seq2))
        aligns = malign.kbest.align(graph, ("0:0", dest), dna_seq1, dna_seq2, 3)

        assert "".join(aligns[0]["a"]) == 'TGGAC-CCGG-G-AAGGTGACCCAC'
        assert "".join(aligns[0]["b"]) == 'TTACCACCGGCGCGAACCCCCCCCC'
        assert math.isclose(aligns[0]["score_a"], 191.0)
        assert math.isclose(aligns[0]["score_b"], 188.0)
        assert math.isclose(aligns[0]["score"], 379.0)

    # TODO: test class only for multiple alignment?
    def test_fill_scorer(self):
        """
        Test scorer filling.
        """

        # TODO: test from zero first
        scorer = malign.utils.fill_matrix("AGCT", "AGCT", malign.utils.DNA_MATRIX)
        assert math.isclose(scorer["A", "A"], malign.utils.DNA_MATRIX["A", "A"])
        assert math.isclose(scorer["A", "T"], malign.utils.DNA_MATRIX["A", "T"])
        assert math.isclose(scorer["A", "-"], malign.utils.DNA_MATRIX["A", "-"])
        assert scorer["-", "-"] == 0.0

        # Remove some entries and check computation from mean values
        del scorer["A", "A"]
        del scorer["A", "T"]
        del scorer["A", "-"]

        scorer = malign.utils.fill_matrix("AGCT", "AGCT", scorer)
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
        graph = malign.kbest.compute_graph(
            dna_seq1, dna_seq2, malign.utils.DNA_MATRIX
        )
        assert len(graph.nodes) == 56
        assert len(graph.edges) == 139
        assert graph.edges["0:0", "1:1"]["weight"] == 2
        assert graph.edges["6:5", "7:6"]["weight"] == 1


if __name__ == "__main__":
    # Explicitly creating and running a test suite allows to profile it
    suite = unittest.TestLoader().loadTestsFromTestCase(TestMalign)
    unittest.TextTestRunner(verbosity=2).run(suite)
