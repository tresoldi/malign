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
        id_matrix = malign.utils.identity_matrix(seqs, match=2, gap_score=-3)

        assert isclose(id_matrix["-", "-", "a"], -3.0)
        assert isclose(id_matrix["a", "a", "a"], 2.0)
        assert isclose(id_matrix["b", "b", "a"], -0.5)

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
            [seq_a, seq_b], k=2, method="anw", matrix=malign.utils.DNA_MATRIX
        )
        assert tuple(alms[0]["seqs"][1]) == ("-", "-", "-", "-", "-", "-", "A")
        assert isclose(alms[0]["score"], -27.0, rel_tol=1e-05)

        seq_a = "GATTACA"
        seq_b = "ATTT"
        alms = malign.multi_align(
            [seq_a, seq_b], k=2, method="anw", matrix=malign.utils.DNA_MATRIX
        )
        assert tuple(alms[0]["seqs"][1]) == ("-", "A", "T", "T", "-", "T", "-")
        assert isclose(alms[0]["score"], 7.0, rel_tol=1e-05)

    def test_nw_alignment_asymmetric(self):
        """
        Test results of alignment with `nw` method with asymmetric matrices.
        """

        # Perform pairwise, assymetric NW
        matrix = malign.utils.DNA_MATRIX.copy()
        matrix["C", "T"] = -99
        seq_a = "GATTACA"
        seq_b = "ATTT"
        alms = malign.multi_align([seq_a, seq_b], k=4, method="anw", matrix=matrix)
        assert tuple(alms[0]["seqs"][1]) == ("-", "A", "T", "T", "T", "-", "-")
        assert isclose(alms[0]["score"], 3.0, rel_tol=1e-05)

    def test_nw_alignment_linguistic(self):
        """
        Test results of alignment with `nw` method on linguistic data.
        """

        filename = Path(__file__).parent.parent
        filename = filename / "docs" / "ita_rus.matrix"
        ita_rus = malign.ScoringMatrix(filename=filename.as_posix())

        alms = malign.multi_align(
            ["Giacomo", "Яков"], k=4, method="anw", matrix=ita_rus
        )
        assert tuple(alms[0]["seqs"][1]) == ("Я", "-", "-", "к", "о", "в", "-")
        assert isclose(alms[0]["score"], 18.0, rel_tol=1e-05)

    def test_multialignment_linguistic(self):
        """
        Test results of alignment with `nw` method on multiwise linguistic data.
        """

        docs_path = Path(__file__).parent.parent
        filename_a = docs_path / "docs" / "ita_rus.matrix"
        filename_b = docs_path / "docs" / "ita_grk.matrix"
        ita_rus = malign.ScoringMatrix(filename=filename_a.as_posix())
        ita_grk = malign.ScoringMatrix(filename=filename_b.as_posix())

        full_matrix = malign.ScoringMatrix(
            scores={}, sub_matrices={(0, 1): ita_rus, (0, 2): ita_grk}
        )
        full_matrix["o", "в", "ο"] = -4
        full_matrix["i", "-", "Ι"] = -4
        full_matrix["c", "к", "κ"] = 10

        seqs = ["Giacomo", "Яков", "Ιακωβος"]
        nw_alms = malign.multi_align(seqs, method="anw", k=4, matrix=full_matrix)
        assert tuple(nw_alms[0]["seqs"][1]) == ("Я", "-", "-", "к", "о", "-", "-", "в")
        assert isclose(nw_alms[0]["score"], 17.8, rel_tol=1e-05)

        yenksp_alms = malign.multi_align(seqs, method="yenksp", k=2, matrix=full_matrix)
        assert tuple(yenksp_alms[0]["seqs"][1]) == ("-", "-", "Я", "к", "о", "в", "-")
        assert isclose(yenksp_alms[0]["score"], 23.458333, rel_tol=1e-05)

    # TODO: reimplement these tests
    def test_alignment_identity(self):
        """
        Test results of alignment with `nw` method on identity matrices.
        """

        seqs = ["VOLDEMORT", "WALDEMAR", "VLADIMIR", "VOLODYMIR"]
        voldemort_matrix = malign.utils.identity_matrix(seqs)

        nw_alms = malign.multi_align(seqs, method="anw", k=4, matrix=voldemort_matrix)
        #        assert tuple(nw_alms[0]["seqs"][0]) == (
        #            "V",
        #            "O",
        #            "L",
        #            "-",
        #            "D",
        #            "E",
        #            "M",
        #            "O",
        #            "R",
        #            "T",
        #        )
        #        assert isclose(nw_alms[0]["score"], -18.654807, rel_tol=1e-05)

        yenksp_alms = malign.multi_align(
            seqs, method="yenksp", k=4, matrix=voldemort_matrix
        )


#        assert tuple(yenksp_alms[0]["seqs"][0]) == (
#            "V",
#            "O",
#            "L",
#            "-",
#            "D",
#            "E",
#            "-",
#            "M",
#            "O",
#            "R",
#            "T",
#        )
#        assert isclose(yenksp_alms[0]["score"], -22.654807, rel_tol=1e-05)


if __name__ == "__main__":
    # Explicitly creating and running a test suite allows to profile it
    suite = unittest.TestLoader().loadTestsFromTestCase(TestMalignResults)
    unittest.TextTestRunner(verbosity=2).run(suite)
