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

# Example cognate sets (DNA sequences that should align)
# Each list contains related sequences we want to align well
cognate_sets = [
    [["A", "C", "G", "T"], ["A", "C", "G", "T"]],  # Perfect match
    [["A", "C", "C", "T"], ["A", "G", "C", "T"]],  # One substitution
    [["T", "G", "A", "C"], ["T", "G", "A", "C"]],  # Perfect match
    [["G", "A", "T", "T"], ["G", "A", "A", "T"]],  # One substitution
]

# Learn optimal matrix using Expectation-Maximization
learned_matrix = malign.learn_matrix(
    cognate_sets,
    method="em",
    max_iter=10,  # Number of EM iterations
    gap="-",
)

# Use the learned matrix for alignment
test_seqs = [["A", "C", "G", "T"], ["A", "G", "G", "T"]]
alignments = malign.align(test_seqs, k=1, matrix=learned_matrix)
print("Learned matrix alignment:")
print(malign.tabulate_alms(alignments))

#' ## Gradient Descent Method
#'
#' Alternative to EM using scipy.optimize for direct optimization:

gd_matrix = malign.learn_matrix(
    cognate_sets,
    method="gradient_descent",
    max_iter=20,  # Optimization iterations
    gap="-",
)

# Compare results
em_alm = malign.align(test_seqs, k=1, matrix=learned_matrix)[0]
gd_alm = malign.align(test_seqs, k=1, matrix=gd_matrix)[0]

print(f"EM score: {em_alm.score:.3f}")
print(f"Gradient descent score: {gd_alm.score:.3f}")

#' ## Using a Custom Initial Matrix
#'
#' You can provide a starting matrix instead of using auto-initialization:

initial_matrix = malign.ScoringMatrix.from_sequences(
    sequences=[["A", "C", "G", "T"], ["A", "C", "G", "T"]],
    match=1.0,
    mismatch=-0.5,
    gap_score=-1.0,
)

refined_matrix = malign.learn_matrix(
    cognate_sets,
    method="em",
    max_iter=5,
    initial_matrix=initial_matrix,
    gap="-",
)

#' ## Phase 2 Notes
#'
#' Current implementation provides basic iteration without convergence checking.
#' The following features will be added in Phase 3:
#' - Convergence criteria and early stopping
#' - Cross-validation for matrix evaluation
#' - Training/validation split utilities
#' - Integration with freqprob for Bayesian priors
#'
#' **Method Selection**:
#' - `method="em"`: Expectation-Maximization, good for interpretability
#' - `method="gradient_descent"`: Direct optimization, potentially faster convergence
#'
#' **Hyperparameters**:
#' - `max_iter`: More iterations allow better convergence (but take longer)
#' - `initial_matrix`: Good initialization can speed up learning
#' - `gap`: Gap symbol must match your sequence encoding
