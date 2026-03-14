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
#' - **Asymmetric scoring**: A->B can have different scores than B->A
#' - **True multi-alignment**: Considers all sequences simultaneously (up to 4)
#' - **k-best results**: Returns multiple high-quality alignments
#' - **Automatic dispatch**: Selects pairwise, N-dim, or progressive alignment
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
#' ## Linguistic Alignment
#'
#' MAlign was designed for computational linguistics. Here's an example
#' comparing cognate words across languages:

# Align phonological segments
sequences = [
    ["k", "a", "t"],  # English "cat"
    ["g", "a", "t", "o"],  # Spanish "gato"
]

alms = malign.align(sequences, k=2, method="anw")
print("\nLinguistic alignment:")
print(malign.tabulate_alms(alms))

#' ## Choosing an Alignment Method
#'
#' MAlign provides two main algorithms:
#'
#' - **`anw`** (default): A* + Needleman-Wunsch. Fast, recommended for most cases.
#' - **`yenksp`**: Yen's k-shortest paths. Maximum quality for small problems.

# Compare methods on the same input
alms_anw = malign.align(["ACGT", "AGCT"], k=3, method="anw")
alms_yen = malign.align(["ACGT", "AGCT"], k=3, method="yenksp")

print("\nANW results:")
print(malign.tabulate_alms(alms_anw))

print("\nYenKSP results:")
print(malign.tabulate_alms(alms_yen))

#' ## K-Best Alignments
#'
#' When the optimal alignment is ambiguous, examining multiple alternatives
#' helps assess uncertainty:

sequences = ["ACGT", "AGT"]
alms = malign.align(sequences, k=5)

print(f"\nTop {len(alms)} alignments:")
for i, aln in enumerate(alms):
    print(f"  #{i}: score={aln.score:.3f}, columns={len(aln.seqs[0])}")

# If top-k scores are similar, the alignment is uncertain
scores = [aln.score for aln in alms]
if len(scores) > 1:
    score_range = max(scores) - min(scores)
    print(
        f"  Score range: {score_range:.3f} {'(uncertain)' if score_range < 0.5 else '(confident)'}"
    )

#' ## Multi-Sequence Alignment
#'
#' MAlign handles 3+ sequences. For up to 4 sequences (small grids),
#' it uses true N-dimensional alignment. For larger sets, it falls
#' back to UPGMA progressive alignment automatically.

sequences = [
    ["k", "a", "t"],
    ["g", "a", "t", "o"],
    ["k", "a", "t", "z", "e"],
]

alms = malign.align(sequences, k=1, method="anw")
print("\n3-sequence alignment:")
print(malign.tabulate_alms(alms))

#' ## Working with Lists of Symbols
#'
#' MAlign works with any list of string symbols -- IPA segments, phonemes,
#' morphemes, etc.:

# IPA segments
seq1 = ["p\u02b0", "a", "t\u02b0"]  # aspirated stops
seq2 = ["b", "a", "d"]

alms = malign.align([seq1, seq2], k=1)
print("\nIPA segment alignment:")
for i, seq in enumerate(alms[0].seqs):
    print(f"  Seq {i}: {list(seq)}")

#' ## Summary
#'
#' - Use `malign.align()` for all alignment tasks
#' - `k` controls how many alternatives to return
#' - `method="anw"` is the default; use `"yenksp"` for maximum quality on small problems
#' - Multi-sequence alignment (3+) is automatic -- no special API needed
#' - `malign.tabulate_alms()` formats results for display
