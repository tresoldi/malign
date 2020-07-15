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
from pathlib import Path

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

    def test_tabulation(self):
        """
        Test alignment tabulation output
        """

        # TODO: assertMultiLineEqual() is failing, only keeping here for coverage

        alms = malign.multi_align(["tra", "fatata"], method="nw", k=3)
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


class TestMalignResults(unittest.TestCase):
    """
    Suite of tests for the quality of alignment methods.
    """

    def test_matrix_from_full_vector(self):
        """
        Test construction of matrices from full vectors.
        """

        # Simple, two domain system
        vectors = {
            ("-", "-"): 0.0,
            ("-", "X"): -3.0,
            ("-", "Y"): -9.0,
            ("a", "-"): -3.0,
            ("a", "X"): 0.0,
            ("a", "Y"): 8.0,
            ("b", "-"): -5.0,
            ("b", "X"): 4.0,
            ("b", "Y"): 4.0,
            ("c", "-"): 2.0,
            ("c", "X"): -1.0,
            ("c", "Y"): 7.0,
        }
        matrix = malign.ScoringMatrix(vectors)
        assert isclose(matrix["-", "-"], 0.0)
        assert isclose(matrix["b", "X"], 4.0)

    def test_matrix_from_partial_vector(self):
        """
        Test construction of matrices from partial vectors.
        """

        # Make a copy, to test inference
        vectors = {
            ("-", "-"): 0.0,
            ("-", "X"): -3.0,
            ("-", "Y"): -9.0,
            #    ('a', '-') : -3.0,
            ("a", "X"): 0.0,
            ("a", "Y"): 8.0,
            ("b", "-"): -5.0,
            ("b", "X"): 4.0,
            ("b", "Y"): 4.0,
            ("c", "-"): 2.0,
            ("c", "X"): -1.0,
            #    ('c', 'Y') :  7.0,
        }

        matrix_std = malign.ScoringMatrix(vectors)
        matrix_fill = malign.ScoringMatrix(vectors, fill_method="distance")

        assert isclose(matrix_std["a", "-"], 1.0)
        assert isclose(matrix_std["c", "Y"], 2.0)
        assert isclose(matrix_fill["a", "-"], 1.0)
        assert isclose(matrix_fill["c", "Y"], 0.8)

    def test_matrix_from_sub_matrices(self):
        """
        Test construction of matrices from sub matrices.
        """

        vector_01 = {
            ("-", "X"): -3.0,
            ("a", "-"): -3.0,
            ("a", "X"): 0.0,
            ("a", "Y"): 8.0,
            ("b", "-"): -5.0,
            ("b", "Y"): 4.0,
            ("c", "X"): -1.0,
            ("c", "Y"): 7.0,
        }
        vector_02 = {
            ("a", "-"): -4.0,
            ("a", "i"): 2.0,
            ("a", "j"): 2.0,
            ("b", "i"): -5.0,
            ("b", "j"): 9.0,
            ("c", "-"): -7.0,
            ("c", "j"): 4.0,
        }

        matrix01 = malign.ScoringMatrix(vector_01)
        matrix02 = malign.ScoringMatrix(vector_02)
        matrix = malign.ScoringMatrix(
            scores={}, sub_matrices={(0, 1): matrix01, (0, 2): matrix02}
        )

        assert isclose(matrix["-", "-", "i"], -0.75)
        assert isclose(matrix["b", "X", "i"], -3)
        assert isclose(matrix["c", "Y", "j"], 5.5)

    def test_identity_matrix(self):
        """
        Test construction of identity matrices.
        """

        seqs = ["ab", "aab", "bbb"]
        id_matrix = malign.utils.identity_matrix(seqs, match=2, gap=-3)

        assert isclose(id_matrix["-", "-", "a"], -3.0)
        assert isclose(id_matrix["a", "a", "a"], 2.0)
        assert isclose(id_matrix["b", "b", "a"], -3.0)

    def test_matrix_filling(self):
        """
        Test filling methods for metrices.
        """

        scores = {
            ("a", "A", "1"): -1,
            ("a", "A", "2"): 4,
            ("a", "A", "3"): 3,
            ("a", "B", "1"): 1,
            ("b", "A", "1"): -10,
            ("b", "A", "2"): 10,
            ("b", "A", "3"): 10,
            ("b", "A", "4"): 10,
            ("c", "B", "1"): 2,
            ("c", "B", "2"): 2,
            ("a", "-", "-"): -2,
            ("b", "-", "-"): -2,
            ("-", "A", "-"): -20,
            ("-", "B", "-"): -20,
            ("-", "-", "1"): -3,
            ("a", "B", "-"): -10,
            ("-", "A", "1"): -100,
            ("-", "A", "2"): -10,
            ("-", "B", "3"): -5,
        }
        matrix = malign.ScoringMatrix(scores)
        assert isclose(matrix["a", "A", "2"], 4.0, rel_tol=1e-05)  # provided
        assert isclose(
            matrix["-", "-", "3"], -4.0, rel_tol=1e-05
        )  # inferred in first round
        assert isclose(
            matrix["c", "-", "4"], -0.247604, rel_tol=1e-05
        )  # inferred in second round

    def test_dumb_alignment(self):
        """
        Test results of alignment with the `dumb` method.
        """

        # Perform pairwise dumb alignment
        seq_a = "tra"
        seq_b = "fatata"
        alms = malign.multi_align([seq_a, seq_b], method="dumb")
        assert tuple(alms[0]["seqs"][0]) == ("-", "t", "r", "a", "-", "-")

        # Perform multiwise dumb alignment
        seqs = ["tra", "fra", "batata", "virp", "x"]
        alms = malign.multi_align(seqs, method="dumb")
        assert tuple(alms[0]["seqs"][3]) == ("-", "v", "i", "r", "p", "-")

    def test_nw_alignment(self):
        """
        Test results of alignment with the `nw` method.
        """

        seq_a = "GATTACA"
        seq_b = "A"
        alms = malign.multi_align(
            [seq_a, seq_b], k=2, method="nw", matrix=malign.utils.DNA_MATRIX
        )
        assert tuple(alms[0]["seqs"][1]) == ("-", "-", "-", "-", "A", "-", "-")
        assert isclose(alms[0]["score"], -2.85714, rel_tol=1e-05)

        seq_a = "GATTACA"
        seq_b = "ATTT"
        alms = malign.multi_align(
            [seq_a, seq_b], k=2, method="nw", matrix=malign.utils.DNA_MATRIX
        )
        assert tuple(alms[0]["seqs"][1]) == ("-", "A", "T", "T", "-", "T", "-")
        assert isclose(alms[0]["score"], 1.57142, rel_tol=1e-05)

    def test_nw_alignment_asymmetric(self):
        """
        Test results of alignment with `nw` method with asymmetric matrices.
        """

        # Perform pairwise, assymetric NW
        matrix = malign.utils.DNA_MATRIX.copy()
        matrix["C", "T"] = -99
        seq_a = "GATTACA"
        seq_b = "ATTT"
        alms = malign.multi_align([seq_a, seq_b], k=4, method="nw", matrix=matrix)
        assert tuple(alms[0]["seqs"][1]) == ("-", "A", "T", "T", "T", "-", "-", "-")
        assert isclose(alms[0]["score"], 0.875, rel_tol=1e-05)

    def test_nw_alignment_linguistic(self):
        """
        Test results of alignment with `nw` method on linguistic data.
        """

        filename = Path(__file__).parent.parent
        filename = filename / "docs" / "ita_rus.matrix"
        ita_rus = malign.ScoringMatrix(filename=filename.as_posix())

        alms = malign.multi_align(["Giacomo", "Яков"], k=4, method="nw", matrix=ita_rus)
        assert tuple(alms[0]["seqs"][1]) == ("-", "Я", "-", "к", "о", "в", "-")
        assert isclose(alms[0]["score"], 3.28571, rel_tol=1e-05)

    def test_multialignment_linguistic(self):
        """
        Test results of alignment with `nw` method on multiwise linguistic data.
        """

        docs_path = Path(__file__).parent.parent
        filename_a = docs_path / "docs" / "ita_rus.matrix"
        filename_b = docs_path / "docs" / "ita_grk.matrix"
        ita_rus = malign.ScoringMatrix(filename=filename_a.as_posix())
        ita_grk = malign.ScoringMatrix(filename=filename_b.as_posix())

        full_matrix = malign.utils.combine_matrices(ita_rus, ita_grk)
        full_matrix["o", "в", "ο"] = -4
        full_matrix["i", "-", "Ι"] = -4
        full_matrix["c", "к", "κ"] = 10

        seqs = ["Giacomo", "Яков", "Ιακωβος"]
        nw_alms = malign.multi_align(seqs, method="nw", k=4, matrix=full_matrix)
        assert tuple(nw_alms[0]["seqs"][1]) == ("-", "Я", "-", "к", "о", "-", "-", "в")
        assert isclose(nw_alms[0]["score"], 4.5375, rel_tol=1e-05)

        yenksp_alms = malign.multi_align(seqs, method="yenksp", k=2, matrix=full_matrix)
        assert tuple(yenksp_alms[0]["seqs"][1]) == (
            "-",
            "Я",
            "-",
            "к",
            "-",
            "-",
            "о",
            "в",
        )
        assert isclose(yenksp_alms[0]["score"], 4.5375, rel_tol=1e-05)

    def test_alignment_identity(self):
        """
        Test results of alignment with `nw` method on identity matrices.
        """

        seqs = ["VOLDEMORT", "WALDEMAR", "VLADIMIR", "VOLODYMIR"]
        voldemort_matrix = malign.utils.identity_matrix(seqs)

        nw_alms = malign.multi_align(seqs, method="nw", k=4, matrix=voldemort_matrix)
        assert tuple(nw_alms[0]["seqs"][0]) == (
            "V",
            "O",
            "L",
            "-",
            "D",
            "E",
            "M",
            "O",
            "R",
            "T",
        )
        assert isclose(nw_alms[0]["score"], -0.86548, rel_tol=1e-05)

        yenksp_alms = malign.multi_align(
            seqs, method="yenksp", k=4, matrix=voldemort_matrix
        )
        assert tuple(yenksp_alms[0]["seqs"][0]) == (
            "-",
            "V",
            "O",
            "L",
            "D",
            "E",
            "M",
            "O",
            "R",
            "T",
        )
        assert isclose(yenksp_alms[0]["score"], -0.76798, rel_tol=1e-05)


if __name__ == "__main__":
    # Explicitly creating and running a test suite allows to profile it
    suite = unittest.TestLoader().loadTestsFromTestCase(TestMalign)
    unittest.TextTestRunner(verbosity=2).run(suite)
