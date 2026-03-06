#' # MAlign Tutorial 2: Scoring Matrices
#'
#' This tutorial explains how to create, use, and customize scoring matrices
#' for alignment operations.
#'
#' ## What are Scoring Matrices?
#'
#' Scoring matrices define the cost/benefit of aligning symbols:
#' - **Match**: High positive score (symbols align well)
#' - **Mismatch**: Lower or negative score (symbols don't match)
#' - **Gap**: Penalty for inserting gaps
#' - **Asymmetric**: A->B can differ from B->A
#'
#' ## Creating Matrices

import malign

#' ### Default (Identity Matrix)
#' If no matrix is provided, MAlign creates an identity matrix automatically:

sequences = ["ACGT", "AGCT"]
alms = malign.align(sequences, k=1)  # Uses default identity matrix
print("Default matrix alignment:")
print(malign.tabulate_alms(alms))

#' ### From Sequences
#'
#' Create a matrix with simple match/mismatch scoring:

matrix = malign.ScoringMatrix.from_sequences(
    sequences=[["A", "C", "G", "T"], ["A", "C", "G", "T"]],
    match=2.0,
    mismatch=-1.0,
    gap_score=-1.5,
)

alms = malign.align(sequences, k=1, matrix=matrix)
print("\nCustom match/mismatch matrix:")
print(malign.tabulate_alms(alms))

#' ### Cross-Domain Matrices
#'
#' When aligning sequences from different alphabets (e.g., Latin and Cyrillic),
#' the two domains have different symbol sets:

matrix = malign.ScoringMatrix.from_sequences(
    sequences=[["A", "C", "G", "T"], ["\u0410", "\u0412", "\u0413", "\u0422"]],
    match=1.0,
    mismatch=-0.5,
    gap_score=-1.0,
)

print(f"\nCross-domain matrix: {matrix.num_domains} domains")
print(f"Domain 0 symbols: {list(matrix.domains[0])}")
print(f"Domain 1 symbols: {list(matrix.domains[1])}")

#' ### From Phonological Features (distfeat)
#'
#' The `from_distfeat()` factory builds matrices from phonological feature
#' distances. Sounds that share more features get higher alignment scores.
#' Requires `pip install malign[features]`.

try:
    matrix = malign.ScoringMatrix.from_distfeat(
        sequences=[["p", "t", "k", "b", "d", "g"], ["p", "t", "k", "b", "d", "g"]],
        gap="-",
        gap_score=-1.0,
    )

    # p-b differ only in voicing: high similarity
    print(f"\nFeature-based scores:")
    print(f"  p-b (voicing only): {matrix[('p', 'b')]:.3f}")
    print(f"  p-d (voicing+place): {matrix[('p', 'd')]:.3f}")
    print(f"  p-g (voicing+place): {matrix[('p', 'g')]:.3f}")
    print(f"  p-p (identical): {matrix[('p', 'p')]:.3f}")

    # Note: from_distfeat() is symmetric
    print(f"  Symmetric: p-b == b-p? {matrix[('p', 'b')] == matrix[('b', 'p')]}")
except ImportError:
    print("\n(Skipping distfeat example -- install with: pip install malign[features])")

#' ### From Substitution Counts
#'
#' When you have observed substitution frequencies (e.g., from a corpus of
#' cognate pairs), use `from_substitution_counts()` to build an asymmetric
#' log-odds matrix:

counts = {
    ("p", "b"): 15,   # p -> b observed 15 times (voicing, common)
    ("b", "p"): 3,    # b -> p observed 3 times (devoicing, rare)
    ("p", "p"): 20,   # p -> p (identity)
    ("b", "b"): 18,   # b -> b (identity)
    ("t", "d"): 10,   # t -> d (voicing)
    ("d", "t"): 2,    # d -> t (devoicing, rare)
    ("t", "t"): 25,
    ("d", "d"): 20,
}
matrix = malign.ScoringMatrix.from_substitution_counts(counts)

print("\nAsymmetric log-odds scores from counts:")
print(f"  p->b (common): {matrix[('p', 'b')]:.3f}")
print(f"  b->p (rare):   {matrix[('b', 'p')]:.3f}")
print(f"  Asymmetric: {matrix[('p', 'b')] != matrix[('b', 'p')]}")

#' ### From YAML Files
#'
#' Matrices can be saved to and loaded from YAML files for reuse:

# Save a matrix
matrix = malign.ScoringMatrix.from_sequences(
    sequences=[["A", "C", "G"], ["A", "C", "G"]],
    match=1.0, mismatch=-0.5, gap_score=-1.0,
)
# matrix.save("my_matrix.yml")

# Load it back
# loaded = malign.ScoringMatrix.from_yaml("my_matrix.yml")

#' ### Manual Construction
#'
#' For complete control, construct a ScoringMatrix directly:

scores = {
    ("A", "A"): 2.0,
    ("A", "C"): -1.0,
    ("A", "-"): -2.0,
    ("C", "A"): -0.5,   # Asymmetric: C->A differs from A->C
    ("C", "C"): 2.0,
    ("C", "-"): -2.0,
    ("-", "A"): -2.0,
    ("-", "C"): -2.0,
    ("-", "-"): 0.0,
}

matrix = malign.ScoringMatrix(
    scores=scores,
    domains=[["-", "A", "C"], ["-", "A", "C"]],
    gap="-",
    impute_method=None,
)

print(f"\nManual matrix: {len(matrix.scores)} score entries")
print(f"  A->C = {matrix[('A', 'C')]:.1f}")
print(f"  C->A = {matrix[('C', 'A')]:.1f} (asymmetric!)")

#' ## Matrix Imputation
#'
#' Sparse matrices (with missing symbol pairs) can be filled automatically:

matrix = malign.ScoringMatrix(
    scores={("A", "A"): 1.0, ("C", "C"): 1.0},
    domains=[["-", "A", "C", "G"], ["-", "A", "C", "G"]],
    gap="-",
    impute_method="mean",  # also "median" or "zero"
)

# Missing scores are imputed on access
print(f"\nImputed matrix:")
print(f"  A-A (known): {matrix[('A', 'A')]:.3f}")
print(f"  A-C (imputed): {matrix[('A', 'C')]:.3f}")

#' ## Inspecting Matrices
#'
#' Use `tabulate()` to view the full matrix:

matrix = malign.ScoringMatrix.from_sequences(
    sequences=[["A", "C", "G"], ["A", "C", "G"]],
    match=1.0, mismatch=-0.5, gap_score=-1.0,
)

print("\nMatrix table:")
print(matrix.tabulate())

#' ## Summary
#'
#' | Factory Method | Symmetric? | Input | Use Case |
#' |---|---|---|---|
#' | `from_sequences()` | Yes | Symbol lists | Quick testing |
#' | `from_distfeat()` | Yes | IPA segments | Phonological knowledge |
#' | `from_substitution_counts()` | No | Frequency counts | Observed sound changes |
#' | `from_yaml()` | Either | YAML file | Reuse saved matrices |
#' | `ScoringMatrix()` | Either | Score dict | Full manual control |
