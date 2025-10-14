#' # MAlign Tutorial 1: Basic Alignment
#'
#' This tutorial introduces the core concepts of MAlign and demonstrates
#' basic sequence alignment operations.
#'
#' ## Installation
#'
#' ```bash
#' pip install malign
#' ```
#'
#' ## Overview
#'
#' MAlign performs multiple sequence alignment with support for:
#' - **Asymmetric scoring**: A→B can have different scores than B→A
#' - **True multi-alignment**: Considers all sequences simultaneously
#' - **k-best results**: Returns multiple high-quality alignments
#'
#' ## Basic DNA Alignment
#'
#' Let's start with a simple example aligning two DNA sequences:

import malign

# Define sequences to align
sequences = ["ATTCGGAT", "TACGGATTT"]

# Perform alignment (returns top 3 alignments)
alignments = malign.align(sequences, k=3, method="anw")

# Display results
print(malign.tabulate_alms(alignments))

#' ## Understanding the Results
#'
#' The output shows:
#' - **Idx**: Alignment index (0 = best, 1 = second best, etc.)
#' - **Seq**: Sequence identifier (A, B, C, ...)
#' - **Score**: Alignment quality score (higher is better)
#' - **Columns**: Aligned positions (- indicates gaps)
#'
#' ## TODO: Complete in Phase 4
#' - Add more examples
#' - Explain scoring
#' - Show different methods (anw vs yenksp)
#' - Demonstrate k-best interpretation
