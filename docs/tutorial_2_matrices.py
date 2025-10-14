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
#' - **Asymmetric**: A→B can differ from B→A
#'
#' ## Creating Matrices

import malign

#' ### Identity Matrix (automatic)
#' If no matrix is provided, MAlign creates an identity matrix:

sequences = ["ACGT", "AGCT"]
alms = malign.align(sequences, k=1)  # Uses default identity matrix

#' ### From Sequences

matrix = malign.ScoringMatrix.from_sequences(
    sequences=[["A", "C", "G", "T"], ["А", "В", "Г", "Т"]],
    match=1.0,
    mismatch=-0.5,
    gap=-1.0
)

#' ### From YAML File (Phase 2)

# matrix = malign.ScoringMatrix.from_yaml("path/to/matrix.yml")

#' ### From Scores (programmatic)

scores = {
    ("A", "А"): 1.0,
    ("A", "В"): -0.5,
    ("A", "Г"): -0.3,
    # ... more scores
}
# matrix = malign.ScoringMatrix.from_scores(scores, domains=[...], gap="-")

#' ## TODO: Complete in Phase 4
#' - Demonstrate matrix imputation
#' - Show asymmetric scoring examples
#' - Explain YAML format
#' - Integration with freqprob/asymcat
