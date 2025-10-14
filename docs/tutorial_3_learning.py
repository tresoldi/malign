#' # MAlign Tutorial 3: Matrix Learning
#'
#' This tutorial demonstrates how to learn optimal scoring matrices
#' from cognate sets or parallel data.
#'
#' ## The Matrix Learning Problem
#'
#' Given cognate sets (sequences we believe are related), we want to find
#' the scoring matrix that maximizes alignment quality across all sets.
#'
#' **Core Assumption**: If cognate sets are correct, their alignment scores
#' should be collectively maximized.
#'
#' ## Basic Matrix Learning (Phase 2 - EM Method)

import malign

# Example cognate sets (Italian-Russian name pairs)
cognate_sets = [
    ["Giacomo", "Яков"],
    ["Pietro", "Пётр"],
    ["Giovanni", "Иван"],
    # ... more cognates
]

# Learn optimal matrix using Expectation-Maximization
# learned_matrix = malign.learn_matrix(
#     cognate_sets,
#     method="em",
#     max_iterations=100,
#     convergence_threshold=0.001
# )

#' ## Cross-Validation (Phase 2)

# Evaluate matrix learning with cross-validation
# cv_results = malign.cross_validate_matrix(
#     cognate_sets,
#     k_folds=5,
#     method="em"
# )
#
# print(f"Mean score: {cv_results['mean_score']:.3f}")
# print(f"Std dev: {cv_results['std_score']:.3f}")

#' ## Integration with External Tools (Phase 2)

# from freqprob import FreqProbDist
#
# # Use frequency distribution as Bayesian prior
# fp_dist = FreqProbDist.from_sequences(all_sequences)
# matrix = malign.learn_matrix(
#     cognate_sets,
#     method="em",
#     prior=fp_dist  # Regularize with frequency prior
# )

#' ## TODO: Complete in Phase 4
#' - Add working examples once learn_matrix() is implemented
#' - Demonstrate EM vs gradient descent
#' - Show convergence plots
#' - Explain hyperparameter tuning
