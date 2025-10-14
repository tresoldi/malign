#' # MAlign Tutorial 4: Advanced Features
#'
#' This tutorial covers advanced features including batch processing,
#' validation metrics, and performance optimization.
#'
#' ## Batch Processing

import malign

# Align multiple cognate sets in batch
# cognate_sets = [
#     ["word1_a", "word1_b"],
#     ["word2_a", "word2_b"],
#     ["word3_a", "word3_b"],
# ]
#
# results = malign.align_batch(
#     cognate_sets,
#     k=5,
#     method="anw",
#     matrix=scorer
# )

#' ## Validation Metrics

# Once Alignment methods are implemented:
# alm = alignments[0]
#
# # Sum-of-pairs score
# sop_score = alm.sum_of_pairs_score(matrix)
#
# # Per-site entropy
# entropy = alm.entropy()
#
# # Consensus sequence
# consensus = alm.consensus(threshold=0.5)
#
# # Gap statistics
# gap_stats = alm.gap_statistics()
# print(f"Total gaps: {gap_stats['total_gaps']}")
# print(f"Gap percentage: {gap_stats['gap_percentage']:.1f}%")

#' ## Algorithm Selection: ANW vs YenKSP
#'
#' ### When to use ANW (default)
#' - Exact alignments needed
#' - Smaller k values (k < 10)
#' - Standard use cases
#'
#' ### When to use YenKSP
#' - Need diverse alignments
#' - Larger k values (k > 10)
#' - Exploring alignment space

# Example with YenKSP
# alms_yenksp = malign.align(sequences, k=20, method="yenksp")

#' ## Performance Considerations
#'
#' **Computational Limits** (to be benchmarked in Phase 3):
#' - Recommended: 2-20 sequences
#' - Maximum practical: ~50 sequences (depends on length)
#' - For 100+ sequences, consider progressive methods
#'
#' ## TODO: Complete in Phase 4
#' - Add performance benchmarks
#' - Demonstrate batch API once implemented
#' - Add validation metrics examples
#' - Include profiling examples
