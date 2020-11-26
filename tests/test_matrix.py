#!/usr/bin/env python3
# pylint: disable=no-self-use

"""
test_malign
===========

Tests for the scoring matrices of the `malign` package.
"""

# TODO: add test for identity matrix
# TODO: add test for initialization only from sparse subdomain
# TODO: add test providing domains
# TODO: add, in general, tests where there is disagreement between scores/subm/domain
# TODO: replace .num_domains with len(.domains) -- or maybe just __len__?

# Import Python libraries
import math
import tempfile
import unittest

# Impor the library itself
import malign

# Vectors for tests
PAIRWISE_TEST_VECTORS = {
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

MULTIWISE_TEST_VECTORS = {
    ("-", "-", "-"): 0.0,
    ("-", "-", "i"): -4.0,
    ("-", "-", "j"): -8.0,
    ("-", "X", "-"): -5.0,
    ("-", "X", "i"): -3.0,
    ("-", "X", "j"): -5.0,
    ("-", "Y", "-"): -5.0,
    ("-", "Y", "i"): -9.0,
    ("-", "Y", "j"): -6.0,
    ("a", "-", "-"): 0.0,
    ("a", "-", "i"): -3.0,
    ("a", "-", "j"): 3.0,
    ("a", "X", "-"): 0.0,
    ("a", "X", "i"): 0.0,
    ("a", "X", "j"): 8.0,
    ("a", "Y", "-"): 8.0,
    ("a", "Y", "i"): 8.0,
    ("a", "Y", "j"): 8.0,
    ("b", "-", "-"): 0.0,
    ("b", "-", "i"): -5.0,
    ("b", "-", "j"): -6.0,
    ("b", "X", "-"): 4.0,
    ("b", "X", "i"): 4.0,
    ("b", "X", "j"): 5.0,
    ("b", "Y", "-"): 4.0,
    ("b", "Y", "i"): 4.0,
    ("b", "Y", "j"): 4.0,
    ("c", "-", "-"): 0.0,
    ("c", "-", "i"): 2.0,
    ("c", "-", "j"): -5.0,
    ("c", "X", "-"): -1.0,
    ("c", "X", "i"): -1.0,
    ("c", "X", "j"): 6.0,
    ("c", "Y", "-"): -7.0,
    ("c", "Y", "i"): 7.0,
    ("c", "Y", "j"): 7.0,
}

PAIRWISE_TEST_SPARSE_VECTOR_01 = {
    ("-", "X"): -3.0,
    ("a", "-"): -3.0,
    ("a", "X"): 0.0,
    ("a", "Y"): 8.0,
    ("b", "-"): -5.0,
    ("b", "Y"): 4.0,
    ("c", "X"): -1.0,
    ("c", "Y"): 7.0,
}
PAIRWISE_TEST_SPARSE_VECTOR_02 = {
    ("a", "-"): -4.0,
    ("a", "i"): 2.0,
    ("a", "j"): 2.0,
    ("b", "i"): -5.0,
    ("b", "j"): 9.0,
    ("c", "-"): -7.0,
    ("c", "j"): 4.0,
}


class TestMalign(unittest.TestCase):
    """
    Suite of tests for scoring matrices.
    """

    def test_pairwise_from_full_vectors(self):
        """
        Test pairwise matrices built from complete vectors.
        """

        # Build matrix
        matrix = malign.ScoringMatrix(PAIRWISE_TEST_VECTORS)

        # Assertions
        assert matrix.num_domains == 2
        assert matrix.gap == "-"
        assert len(matrix.scores) == 12
        assert matrix["-", "-"] == 0.0
        assert matrix["a", "Y"] == 8.0
        assert len(matrix.domains) == 2
        assert tuple(matrix.domains[1]) == ("-", "X", "Y")

    def test_pairwise_from_full_vectors_with_domains(self):
        """
        Test pairwise matrices built from complete vectors with domains.
        """

        # Build matrix with "correct" domains
        matrix_a = malign.ScoringMatrix(
            PAIRWISE_TEST_VECTORS, domains=[["-", "a", "b", "c"], ["-", "X", "Y"]]
        )

        # Build matrix with "expanded" domains
        matrix_b = malign.ScoringMatrix(
            PAIRWISE_TEST_VECTORS,
            domains=[["-", "a", "b", "c", "d"], ["-", "X", "Y", "Z"]],
        )

        # Build matrix with "insufficient" domains
        with self.assertRaises(ValueError):
            malign.ScoringMatrix(
                PAIRWISE_TEST_VECTORS, domains=[["-", "a", "b"], ["-", "Y", "Z"]]
            )

        # Assertions
        assert tuple(matrix_a.domains[1]) == ("-", "X", "Y")
        assert tuple(matrix_b.domains[1]) == ("-", "X", "Y", "Z")

    def test_multiwise_from_full_vectors(self):
        """
        Test multiwise matrices built from complete vectors.
        """

        # Build matrix
        matrix = malign.ScoringMatrix(MULTIWISE_TEST_VECTORS)

        # Assertions
        assert matrix.num_domains == 3
        assert matrix.gap == "-"
        assert len(matrix.scores) == 69
        assert len(matrix.domains) == 3
        assert tuple(matrix.domains[2]) == ("-", "i", "j")
        assert matrix["-", "-", "-"] == 0.0
        assert math.isclose(matrix["a", "Y", "j"], 8.0)

    def test_multiwise_from_sparse_vectors(self):
        """
        Test multiwise matrices built from sparse vectors.
        """

        # Build matrices with the various filling methods
        vectors = MULTIWISE_TEST_VECTORS.copy()
        vectors.pop(("a", "X", "j"))
        vectors.pop(("c", "X", "-"))
        vectors.pop(("b", "-", "i"))
        vectors.pop(("b", "-", "j"))
        vectors.pop(("c", "Y", "j"))
        vectors.pop(("-", "X", "-"))
        matrix_default = malign.ScoringMatrix(vectors)
        matrix_dt = malign.ScoringMatrix(vectors, impute_method="decision_tree")
        matrix_et = malign.ScoringMatrix(vectors, impute_method="extra_trees")
        matrix_kn = malign.ScoringMatrix(vectors, impute_method="k_neighbors")
        matrix_br = malign.ScoringMatrix(vectors, impute_method="bayesian_ridge")

        # Assertions
        assert matrix_default.num_domains == 3
        assert matrix_default.gap == "-"
        assert len(matrix_default.scores) == 69
        assert len(matrix_default.domains) == 3
        assert tuple(matrix_default.domains[2]) == ("-", "i", "j")

        assert matrix_default["-", "-", "-"] == 0.0
        assert math.isclose(matrix_default["a", "Y", "j"], 8.0, rel_tol=1e-03)
        assert math.isclose(matrix_default["a", "X", "j"], 2.83974, rel_tol=1e-03)
        assert math.isclose(matrix_default["c", "X", "-"], 0.48655, rel_tol=1e-03)
        assert math.isclose(matrix_default["b", "-", "i"], 1.36074, rel_tol=1e-03)
        assert math.isclose(matrix_default["b", "-", "j"], 1.45179, rel_tol=1e-03)
        assert math.isclose(matrix_default["c", "Y", "j"], 1.43649, rel_tol=1e-03)
        assert math.isclose(matrix_default["-", "X", "-"], -3.45639, rel_tol=1e-03)

        assert math.isclose(matrix_dt["a", "Y", "j"], 8.0, rel_tol=1e-03)
        assert math.isclose(matrix_dt["a", "X", "j"], 0.0, rel_tol=1e-03)
        assert math.isclose(matrix_dt["c", "X", "-"], 0.0, rel_tol=1e-03)
        assert math.isclose(matrix_dt["b", "-", "i"], -4.0, rel_tol=1e-03)
        assert math.isclose(matrix_dt["b", "-", "j"], -8.0, rel_tol=1e-03)
        assert math.isclose(matrix_dt["c", "Y", "j"], -6.0, rel_tol=1e-03)
        assert math.isclose(matrix_dt["-", "X", "-"], 0.0, rel_tol=1e-03)

        assert math.isclose(matrix_et["a", "Y", "j"], 8.0, rel_tol=1e-03)
        assert math.isclose(matrix_et["a", "X", "j"], 0.0, rel_tol=1e-03)
        assert math.isclose(matrix_et["c", "X", "-"], -7.0, rel_tol=1e-03)
        assert math.isclose(matrix_et["b", "-", "i"], -3.0, rel_tol=1e-03)
        assert math.isclose(matrix_et["b", "-", "j"], 3.0, rel_tol=1e-03)
        assert math.isclose(matrix_et["c", "Y", "j"], 7.0, rel_tol=1e-03)
        assert math.isclose(matrix_et["-", "X", "-"], -4.5, rel_tol=1e-03)

        assert math.isclose(matrix_kn["a", "Y", "j"], 8.0, rel_tol=1e-03)
        assert math.isclose(matrix_kn["a", "X", "j"], 1.86666, rel_tol=1e-03)
        assert math.isclose(matrix_kn["c", "X", "-"], 0.93333, rel_tol=1e-03)
        assert math.isclose(matrix_kn["b", "-", "i"], 1.46666, rel_tol=1e-03)
        assert math.isclose(matrix_kn["b", "-", "j"], 1.2, rel_tol=1e-03)
        assert math.isclose(matrix_kn["c", "Y", "j"], 1.6, rel_tol=1e-03)
        assert math.isclose(matrix_kn["-", "X", "-"], 0.06666, rel_tol=1e-03)

        assert math.isclose(matrix_br["a", "Y", "j"], 8.0, rel_tol=1e-03)
        assert math.isclose(matrix_br["a", "X", "j"], 2.83974, rel_tol=1e-03)
        assert math.isclose(matrix_br["c", "X", "-"], 0.48655, rel_tol=1e-03)
        assert math.isclose(matrix_br["b", "-", "i"], 1.36074, rel_tol=1e-03)
        assert math.isclose(matrix_br["b", "-", "j"], 1.45179, rel_tol=1e-03)
        assert math.isclose(matrix_br["c", "Y", "j"], 1.43649, rel_tol=1e-03)
        assert math.isclose(matrix_br["-", "X", "-"], -3.45639, rel_tol=1e-03)

    def test_multiwise_from_subvectors(self):
        """
        Test multiwise matrices built from sub vectors.
        """

        # Build sub matrices, and then the main matrix
        scores_01 = {
            (key[0], key[1], None): value
            for key, value in PAIRWISE_TEST_SPARSE_VECTOR_01.items()
        }
        scores_02 = {
            (key[0], None, key[1]): value
            for key, value in PAIRWISE_TEST_SPARSE_VECTOR_02.items()
        }
        scores = {**scores_01, **scores_02}
        matrix = malign.ScoringMatrix(scores)

        # Assertions
        assert matrix.num_domains == 3
        assert matrix.gap == "-"
        assert len(matrix.scores) == 69
        assert len(matrix.domains) == 3

        assert matrix["-", "-", "-"] == 0.0
        assert math.isclose(matrix["a", "Y", "j"], 7.29478, rel_tol=1e-05)
        assert math.isclose(matrix["a", "X", "j"], 2.91425, rel_tol=1e-05)
        assert math.isclose(matrix["c", "X", "-"], -2.40041, rel_tol=1e-05)
        assert math.isclose(matrix["b", "-", "i"], -1.67591, rel_tol=1e-05)
        assert math.isclose(matrix["b", "-", "j"], 1.87231, rel_tol=1e-05)
        assert math.isclose(matrix["c", "Y", "j"], 6.80505, rel_tol=1e-05)
        assert math.isclose(matrix["-", "X", "-"], -2.19299, rel_tol=1e-05)

    def test_subdomain_query(self):
        """
        Test querying of subdomains.
        """

        # Build matrices with the various filling methods
        matrix = malign.ScoringMatrix(MULTIWISE_TEST_VECTORS)

        assert math.isclose(matrix[None, "X", "i"], 0.653183, rel_tol=1e-05)
        assert math.isclose(matrix["c", None, "i"], 0.53235, rel_tol=1e-05)
        assert math.isclose(matrix["c", "X", None], 1.33871, rel_tol=1e-05)

    def test_load_save(self):
        """
        Test load and saving matrices
        """

        # Build matrices with the various filling methods
        matrix = malign.ScoringMatrix(MULTIWISE_TEST_VECTORS)

        # Build a temporary file name and save
        handler = tempfile.NamedTemporaryFile()
        matrix.save(handler.name)

        # Load and check
        matrix2 = malign.ScoringMatrix(handler.name)

        # Assertions
        assert matrix.scores == matrix2.scores
        assert tuple(matrix.domains) == tuple(matrix2.domains)

    def test_copy(self):
        """
        Test method for matrix copy.
        """

        # Build reference matrix
        ref_matrix = malign.ScoringMatrix(MULTIWISE_TEST_VECTORS)

        # Get copy
        cpy_matrix = ref_matrix.copy()

        # Perform manual comparison
        assert ref_matrix.scores == cpy_matrix.scores
        assert ref_matrix.domains == cpy_matrix.domains

    def test_set_item(self):
        """
        Test matrix __setitem__.
        """

        # Build reference matrix
        matrix = malign.ScoringMatrix(MULTIWISE_TEST_VECTORS)

        # Various sets and tests
        matrix["a", "X", "i"] = -111
        matrix[None, "X", "i"] = -222
        with self.assertRaises(ValueError):
            matrix["<", "X", "i"] = -333

        assert matrix["a", "X", "i"] == -111
        assert matrix[None, "X", "i"] == -222

    def test_tabulate(self):
        """
        Test matrix tabulation.
        """

        # Build reference matrix
        matrix_a = malign.ScoringMatrix(PAIRWISE_TEST_VECTORS)
        matrix_b = malign.ScoringMatrix(MULTIWISE_TEST_VECTORS)

        # NOTE: currently only building it, to get coverage
        assert len(matrix_a.tabulate()) > 0
        assert len(matrix_b.tabulate()) > 0


if __name__ == "__main__":
    # Explicitly creating and running a test suite allows to profile it
    suite = unittest.TestLoader().loadTestsFromTestCase(TestMalign)
    unittest.TextTestRunner(verbosity=2).run(suite)
