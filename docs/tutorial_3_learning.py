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
#' MAlign offers three learning approaches:
#' - **`learn_matrix()`**: Supervised, requires pre-clustered cognate sets
#' - **`bootstrap_matrix()`**: Unsupervised, works with arbitrary pairs
#' - **Prior-guided bootstrap**: Blend phonological priors with data
#'
#' ## Supervised Learning: learn_matrix()
#'
#' ### EM Method

import malign

# Example cognate sets -- each list contains related sequences
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
    max_iter=10,
    gap="-",
    verbose=True,
)

# Use the learned matrix for alignment
test_seqs = [["A", "C", "G", "T"], ["A", "G", "G", "T"]]
alignments = malign.align(test_seqs, k=1, matrix=learned_matrix)
print("\nLearned matrix alignment:")
print(malign.tabulate_alms(alignments))

#' ### Gradient Descent Method
#'
#' More accurate with limited data, uses L-BFGS-B optimization:

gd_matrix = malign.learn_matrix(
    cognate_sets,
    method="gradient_descent",
    max_iter=20,
    gap="-",
    bounds=(-10.0, 10.0),
    patience=5,
)

# Compare results
em_alm = malign.align(test_seqs, k=1, matrix=learned_matrix)[0]
gd_alm = malign.align(test_seqs, k=1, matrix=gd_matrix)[0]

print(f"\nEM score: {em_alm.score:.3f}")
print(f"Gradient descent score: {gd_alm.score:.3f}")

#' ### Convergence Features
#'
#' Both methods support:
#' - **Score-based convergence**: Stops when relative score change < threshold
#' - **Matrix-based convergence**: Stops when Frobenius norm change < threshold
#' - **Early stopping**: Patience-based -- stops after N iterations without improvement
#'
#' ### Custom Initial Matrix
#'
#' Providing a good starting matrix can improve convergence:

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

#' ## Unsupervised Learning: bootstrap_matrix()
#'
#' When you have sequence pairs but no pre-clustered cognate sets,
#' use `bootstrap_matrix()`. It iteratively aligns pairs and
#' re-estimates scores using log-odds with Laplace smoothing.

# Simulated phonological pairs: voicing changes (p->b, t->d, k->g)
pairs = [
    (["p", "a", "t", "a"], ["b", "a", "d", "a"]),
    (["p", "a", "k", "a"], ["b", "a", "g", "a"]),
    (["t", "a", "p", "a"], ["d", "a", "b", "a"]),
    (["k", "a", "t", "a"], ["g", "a", "d", "a"]),
    (["p", "i", "t", "i"], ["b", "i", "d", "i"]),
]

matrix = malign.bootstrap_matrix(pairs, max_iter=10, verbose=True)

print(f"\nBootstrap matrix: {matrix.num_domains} domains, {len(matrix.scores)} scores")

# Check that voicing correspondences score higher than random
print(f"  p->b score: {matrix.scores.get(('p', 'b'), 0):.3f}")
print(f"  p->g score: {matrix.scores.get(('p', 'g'), 0):.3f}")

# Use the learned matrix
alms = malign.align(
    [["p", "a", "t"], ["b", "a", "d"]],
    k=1,
    matrix=matrix,
)
print("\nBootstrap alignment:")
print(malign.tabulate_alms(alms))

#' ## Prior-Guided Learning
#'
#' When you have phonological knowledge (e.g., from `distfeat`) but want
#' to refine it with language-specific data, use `prior_matrix`.
#' The prior weight decays linearly to zero, so phonological structure
#' dominates early and fades as data-driven scores stabilize.

try:
    # Step 1: Build a phonological prior
    prior = malign.ScoringMatrix.from_distfeat(
        sequences=[
            ["p", "t", "k", "b", "d", "g", "a", "i"],
            ["p", "t", "k", "b", "d", "g", "a", "i"],
        ],
    )

    # Step 2: Refine with cognate pairs
    matrix_with_prior = malign.bootstrap_matrix(
        pairs,
        max_iter=10,
        prior_matrix=prior,   # phonological feature prior
        prior_weight=0.5,     # initial blending strength (decays to 0)
    )

    # Step 3: Compare with vs without prior
    matrix_no_prior = malign.bootstrap_matrix(pairs, max_iter=10)

    common_keys = set(matrix_no_prior.scores) & set(matrix_with_prior.scores)
    diffs = [
        abs(matrix_no_prior.scores[k] - matrix_with_prior.scores[k])
        for k in common_keys
    ]
    print(f"\nPrior effect: max score diff = {max(diffs):.3f}")

except ImportError:
    print("\n(Skipping prior example -- install with: pip install malign[features])")

#' ### Tuning prior_weight
#'
#' - `0.3-0.5`: Moderate regularization (recommended)
#' - `0.7-1.0`: Strong regularization (useful with very few pairs)
#' - The weight decays linearly: at iteration `i`, effective weight =
#'   `prior_weight * (1 - i / max_iter)`
#'
#' ## Block-Aware Bootstrap
#'
#' When your data includes diphthongization or metathesis patterns,
#' complementary-gap columns inflate gap-symbol counts. Use `block_merge=True`
#' to reduce this effect:

pairs_with_blocks = [
    (["a", "t", "a"], ["j", "e", "d", "a"]),
    (["p", "a", "t"], ["b", "a", "d"]),
    (["p", "a", "t"], ["b", "a", "d"]),
]

matrix_blocks = malign.bootstrap_matrix(
    pairs_with_blocks,
    max_iter=5,
    block_merge=True,
    max_block_size=2,
)

print(f"\nBlock-aware matrix: {len(matrix_blocks.scores)} scores")

#' ## Method Selection Guide
#'
#' | Situation | Method |
#' |---|---|
#' | Pre-clustered cognate sets, 100+ sets | `learn_matrix(method="em")` |
#' | Pre-clustered cognate sets, few sets | `learn_matrix(method="gradient_descent")` |
#' | Arbitrary sequence pairs | `bootstrap_matrix()` |
#' | Pairs + phonological knowledge | `bootstrap_matrix(prior_matrix=...)` |
#' | Pairs with diphthongization/metathesis | `bootstrap_matrix(block_merge=True)` |
#'
#' ## Saving and Reusing Learned Matrices

matrix = malign.bootstrap_matrix(pairs, max_iter=5)

# Save for later use
# matrix.save("learned_voicing.yml")

# Load and use
# matrix = malign.ScoringMatrix.from_yaml("learned_voicing.yml")
# alms = malign.align(new_sequences, k=3, matrix=matrix)
