"""Tests for matrix learning and optimization (Phase 2)."""


# All tests marked as TODO for Phase 2 when learn_matrix() is implemented

# TODO: Phase 2
# def test_em_learning_basic():
#     """Test basic EM matrix learning."""
#     cognate_sets = [["word1a", "word1b"], ["word2a", "word2b"]]
#     matrix = malign.learn_matrix(cognate_sets, method="em")
#     assert matrix is not None
#
#
# def test_cross_validation():
#     """Test cross-validation of matrix learning."""
#     cognate_sets = [["w1a", "w1b"], ["w2a", "w2b"], ["w3a", "w3b"]]
#     cv_results = malign.cross_validate_matrix(cognate_sets, k_folds=3)
#     assert "mean_score" in cv_results
#     assert "std_score" in cv_results
#
#
# def test_gradient_descent_learning():
#     """Test gradient descent matrix learning."""
#     cognate_sets = [["word1a", "word1b"], ["word2a", "word2b"]]
#     matrix = malign.learn_matrix(cognate_sets, method="gradient_descent")
#     assert matrix is not None
