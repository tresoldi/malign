"""
test_matrix
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
import pytest

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

MULTIWISE_TEST_VECTORS_SPARSE = {
    ("-", "-", "-"): 0.0,
    ("-", "-", "i"): -4.0,
    ("-", "-", "j"): -8.0,
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
    ("a", "Y", "-"): 8.0,
    ("a", "Y", "i"): 8.0,
    ("a", "Y", "j"): 8.0,
    ("b", "-", "-"): 0.0,
    ("b", "X", "-"): 4.0,
    ("b", "X", "i"): 4.0,
    ("b", "X", "j"): 5.0,
    ("b", "Y", "-"): 4.0,
    ("b", "Y", "i"): 4.0,
    ("b", "Y", "j"): 4.0,
    ("c", "-", "-"): 0.0,
    ("c", "-", "i"): 2.0,
    ("c", "-", "j"): -5.0,
    ("c", "X", "i"): -1.0,
    ("c", "X", "j"): 6.0,
    ("c", "Y", "-"): -7.0,
    ("c", "Y", "i"): 7.0,
}


def test_pairwise_from_full_vectors():
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


def test_pairwise_from_full_vectors_with_domains():
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
    with pytest.raises(ValueError):
        malign.ScoringMatrix(
            PAIRWISE_TEST_VECTORS, domains=[["-", "a", "b"], ["-", "Y", "Z"]]
        )

    # Assertions
    assert tuple(matrix_a.domains[1]) == ("-", "X", "Y")
    assert tuple(matrix_b.domains[1]) == ("-", "X", "Y", "Z")


def test_multiwise_from_full_vectors():
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
    assert matrix["a", "Y", "j"] == pytest.approx(8.0)


# TODO: add test with "default" (currently "mean")
@pytest.mark.parametrize(
    "method,num_domains,gap,size,tests",
    [
        [
            "mean",
            3,
            "-",
            69,
            [
                [("a", "Y", "j"), 8.0, 1e-03],
                [("a", "X", "j"), 0.366666, 1e-03],
                [("c", "X", "-"), 0.366666, 1e-03],
                [("b", "-", "i"), 0.366666, 1e-03],
                [("b", "-", "j"), 0.366666, 1e-03],
                [("c", "Y", "j"), 0.366666, 1e-03],
                [("-", "X", "-"), 0.366666, 1e-03],
            ],
        ],
        [
            "median",
            3,
            "-",
            69,
            [
                [("a", "Y", "j"), 8.0, 1e-03],
                [("a", "X", "j"), 0.0, 1e-03],
                [("c", "X", "-"), 0.0, 1e-03],
                [("b", "-", "i"), 0.0, 1e-03],
                [("b", "-", "j"), 0.0, 1e-03],
                [("c", "Y", "j"), 0.0, 1e-03],
                [("-", "X", "-"), 0.0, 1e-03],
            ],
        ],
        [
            "decision_tree",
            3,
            "-",
            69,
            [
                [("a", "Y", "j"), 8.0, 1e-03],
                [("a", "X", "j"), 0.0, 1e-03],
                [("c", "X", "-"), 0.0, 1e-03],
                [("b", "-", "i"), -4.0, 1e-03],
                [("b", "-", "j"), -8.0, 1e-03],
                [("c", "Y", "j"), -6.0, 1e-03],
                [("-", "X", "-"), 0.0, 1e-03],
            ],
        ],
        [
            "extra_trees",
            3,
            "-",
            69,
            [
                [("a", "Y", "j"), 8.0, 1e-03],
                [("a", "X", "j"), 0.0, 1e-03],
                [("c", "X", "-"), -7.0, 1e-03],
                [("b", "-", "i"), -3.0, 1e-03],
                [("b", "-", "j"), 3.0, 1e-03],
                [("c", "Y", "j"), 7.0, 1e-03],
                [("-", "X", "-"), -4.5, 1e-03],
            ],
        ],
        [
            "k_neighbors",
            3,
            "-",
            69,
            [
                [("a", "Y", "j"), 8.0, 1e-03],
                [("a", "X", "j"), 1.86666, 1e-03],
                [("c", "X", "-"), 0.93333, 1e-03],
                [("b", "-", "i"), 1.46666, 1e-03],
                [("b", "-", "j"), 1.2, 1e-03],
                [("c", "Y", "j"), 1.6, 1e-03],
                [("-", "X", "-"), 0.06666, 1e-03],
            ],
        ],
        [
            "bayesian_ridge",
            3,
            "-",
            69,
            [
                [("a", "Y", "j"), 8.0, 1e-03],
                [("a", "X", "j"), 2.83974, 1e-03],
                [("c", "X", "-"), 0.48655, 1e-03],
                [("b", "-", "i"), 1.36074, 1e-03],
                [("b", "-", "j"), 1.45179, 1e-03],
                [("c", "Y", "j"), 1.43649, 1e-03],
                [("-", "X", "-"), -3.45639, 1e-03],
            ],
        ],
    ],
)
def test_multiwise_from_sparse_vectors(method, num_domains, gap, size, tests):
    """
    Test multiwise matrices built from sparse vectors.
    """

    matrix = malign.ScoringMatrix(MULTIWISE_TEST_VECTORS_SPARSE, impute_method=method)

    assert matrix.num_domains == num_domains
    assert matrix.gap == gap
    assert len(matrix.scores) == size
    assert ["-", "i", "j"] in matrix.domains

    assert matrix["-", "-", "-"] == 0.0
    for key, expected, rel in tests:
        assert matrix[key] == pytest.approx(expected, rel=rel)


def test_multiwise_from_subvectors():
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
    assert math.isclose(matrix["a", "Y", "j"], 0.5, rel_tol=1e-05)
    assert math.isclose(matrix["a", "X", "j"], 0.5, rel_tol=1e-05)
    assert math.isclose(matrix["c", "X", "-"], 0.5, rel_tol=1e-05)
    assert math.isclose(matrix["b", "-", "i"], 0.5, rel_tol=1e-05)
    assert math.isclose(matrix["b", "-", "j"], 0.5, rel_tol=1e-05)
    assert math.isclose(matrix["c", "Y", "j"], 0.5, rel_tol=1e-05)
    assert math.isclose(matrix["-", "X", "-"], 0.5, rel_tol=1e-05)


def test_subdomain_query():
    """
    Test querying of subdomains.
    """

    # Build matrices with the various filling methods
    matrix = malign.ScoringMatrix(MULTIWISE_TEST_VECTORS)
    assert math.isclose(matrix[None, "X", "i"], 0.25, rel_tol=1e-05)
    assert math.isclose(matrix["c", None, "i"], 0.25, rel_tol=1e-05)
    assert math.isclose(matrix["c", "X", None], 0.25, rel_tol=1e-05)


def test_load_save():
    """
    Test load and saving matrices
    """

    # Build matrices with the various filling methods
    matrix = malign.ScoringMatrix(MULTIWISE_TEST_VECTORS)

    # Build a temporary file name and save
    # TODO: does not work on windows...

    #        handler = tempfile.NamedTemporaryFile()
    #        matrix.save(handler.name)
    #
    #        # Load and check
    #        matrix2 = malign.ScoringMatrix(handler.name)
    #
    #        # Assertions
    #        assert matrix.scores == matrix2.scores
    #        assert tuple(matrix.domains) == tuple(matrix2.domains)


def test_copy():
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

    # Assert they are different
    assert id(ref_matrix) != id(cpy_matrix)


def test_set_item():
    """
    Test matrix __setitem__.
    """

    # Build reference matrix
    matrix = malign.ScoringMatrix(MULTIWISE_TEST_VECTORS)

    # Various sets and tests
    matrix["a", "X", "i"] = -111
    matrix[None, "X", "i"] = -222
    with pytest.raises(ValueError):
        matrix["<", "X", "i"] = -333

    assert matrix["a", "X", "i"] == -111
    assert matrix[None, "X", "i"] == -222


def test_tabulate():
    """
    Test matrix tabulation.
    """

    # Build reference matrix
    matrix_a = malign.ScoringMatrix(PAIRWISE_TEST_VECTORS)
    matrix_b = malign.ScoringMatrix(MULTIWISE_TEST_VECTORS)

    # NOTE: currently only building it, to get coverage
    assert len(matrix_a.tabulate()) > 0
    assert len(matrix_b.tabulate()) > 0
