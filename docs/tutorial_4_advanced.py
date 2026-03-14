#' # MAlign Tutorial 4: Advanced Features
#'
#' This tutorial covers advanced features including block detection,
#' evaluation metrics, multi-sequence dispatch, and performance tips.
#'
#' ## Block Detection and Column Merging
#'
#' In phonological alignment, certain sound changes produce 1-to-many
#' correspondences that appear as complementary gaps:
#'
#' ```
#' Diphthongization (a -> je):     Metathesis (tr -> rt):
#'   a  -                           t  r  -
#'   j  e                           -  r  t
#' ```
#'
#' MAlign detects these **block patterns** and merges them into compound
#' symbols (tuples).
#'
#' ### Automatic Merging via align()

import malign

# Diphthongization: "a" becomes "je"
alms = malign.align(
    [["a"], ["j", "e"]],
    k=1,
    merge_blocks=True,
)

print("Diphthongization (merge_blocks=True):")
for i, seq in enumerate(alms[0].seqs):
    print(f"  Seq {i}: {list(seq)}")
# Seq 0: ["a"]
# Seq 1: [("j", "e")]  -- compound tuple

#' ### Manual Merging

# Create an alignment with complementary gaps
aln = malign.Alignment(
    seqs=[("a", "-", "x"), ("j", "e", "x")],
    score=1.0,
)

merged = malign.merge_alignment_blocks(aln)
print("\nManual merge:")
print(f"  Before: {[list(s) for s in aln.seqs]}")
print(f"  After:  {[list(s) for s in merged.seqs]}")
print(f"  Score preserved: {merged.score == aln.score}")

#' ### How Block Detection Works
#'
#' 1. Columns classified as FULL (no gaps), PARTIAL (some gaps), ALL_GAP
#' 2. Maximal runs of PARTIAL columns found
#' 3. Runs extended to adjacent FULL columns if within `max_block_size`
#' 4. Within a block per sequence:
#'    - 0 non-gap symbols -> gap symbol
#'    - 1 symbol -> unwrapped
#'    - 2+ symbols -> tuple

from malign.blocks import detect_blocks

# Detect blocks directly
seqs = [("a", "-", "x", "b", "-"), ("j", "e", "x", "-", "c")]
blocks = detect_blocks(seqs)
print(f"\nDetected blocks in 5-column alignment: {blocks}")
# [(0, 1), (3, 4)] -- two 2-column blocks

#' ### Configuring max_block_size
#'
#' The default `max_block_size=2` handles diphthongization and simple
#' metathesis. Increase for longer patterns:

# F P F pattern with max_block_size=2 vs 3
seqs = [("a", "-", "b"), ("j", "e", "b")]
print(f"\nF P F with max_block_size=2: {detect_blocks(seqs, max_block_size=2)}")
print(f"F P F with max_block_size=3: {detect_blocks(seqs, max_block_size=3)}")

#' ## Evaluation Metrics
#'
#' MAlign provides standard alignment evaluation metrics:

# Create a "gold standard" alignment
gold = malign.Alignment(
    seqs=[("A", "C", "G", "T"), ("A", "C", "G", "T")],
    score=4.0,
)

# Create a predicted alignment (one column differs)
predicted = malign.Alignment(
    seqs=[("A", "C", "G", "T"), ("A", "G", "C", "T")],
    score=2.0,
)

accuracy = malign.alignment_accuracy(predicted, gold)
precision, recall = malign.alignment_precision_recall(predicted, gold)
f1 = malign.alignment_f1(predicted, gold)

print("\nEvaluation metrics:")
print(f"  Accuracy:  {accuracy:.2%}")
print(f"  Precision: {precision:.2%}")
print(f"  Recall:    {recall:.2%}")
print(f"  F1:        {f1:.2%}")

#' ## Multi-Sequence Alignment Dispatch
#'
#' MAlign automatically selects the best alignment strategy:
#'
#' | Sequences | Strategy | Quality |
#' |---|---|---|
#' | N=2 | Direct pairwise | Optimal |
#' | N=3-4 (small grid) | True N-dimensional | Globally optimal |
#' | N=5+ or large grid | UPGMA progressive | Near-optimal |
#'
#' The dispatch is transparent -- just call `align()`:

# 2 sequences: direct pairwise
alms2 = malign.align(["ACGT", "AGCT"], k=1)
print(f"\n2-seq alignment: {len(alms2[0].seqs[0])} columns")

# 3 sequences: true N-dimensional (if grid is small enough)
alms3 = malign.align(
    [["k", "a", "t"], ["g", "a", "t", "o"], ["k", "a", "t", "z"]],
    k=1,
)
print(f"3-seq alignment: {len(alms3[0].seqs[0])} columns")

#' ## Algorithm Selection: ANW vs YenKSP
#'
#' ```
#' How many sequences?
#' |-- 2-4 sequences
#' |   |-- k <= 10? -> Use YenKSP (maximum quality)
#' |   |-- k > 10?  -> Use ANW (faster for large k)
#' |-- 5+ sequences -> Use ANW (YenKSP too slow)
#' ```

# Small problem, want maximum quality
alms = malign.align(["cat", "bat"], k=5, method="yenksp")
print(f"\nYenKSP (2 seqs, k=5): {len(alms)} alignments")

# Large k: use ANW
alms = malign.align(["cat", "bat"], k=20, method="anw")
print(f"ANW (2 seqs, k=20): {len(alms)} alignments")

#' ## Full Pipeline Example
#'
#' Combining multiple features for a realistic linguistic workflow:

# Step 1: Learn a matrix from cognate pairs
pairs = [
    (["p", "a", "t", "a"], ["b", "a", "d", "a"]),
    (["p", "a", "k", "a"], ["b", "a", "g", "a"]),
    (["t", "a", "p", "a"], ["d", "a", "b", "a"]),
    (["k", "a", "t", "a"], ["g", "a", "d", "a"]),
]

matrix = malign.bootstrap_matrix(pairs, max_iter=10)

# Step 2: Align new data with learned matrix
alms = malign.align(
    [["p", "a", "t"], ["b", "a", "d"]],
    k=3,
    matrix=matrix,
    merge_blocks=True,
)

print("\nFull pipeline:")
print(f"  Learned matrix with {len(matrix.scores)} scores")
print(f"  Got {len(alms)} alignments")
print(malign.tabulate_alms(alms[:2]))

#' ## Performance Tips
#'
#' 1. **Choose the right algorithm**: Use `anw` for 5+ sequences or k > 10
#' 2. **Limit k**: Don't compute more alternatives than needed
#' 3. **Use learned matrices**: They often produce faster convergence
#' 4. **Batch with multiprocessing**: Align multiple sets in parallel
#'
#' ```python
#' from multiprocessing import Pool
#'
#' def align_task(sequences):
#'     return malign.align(sequences, k=3, method="anw")
#'
#' with Pool(4) as pool:
#'     results = pool.map(align_task, sequence_batches)
#' ```
#'
#' ## Summary
#'
#' | Feature | API |
#' |---|---|
#' | Block merging | `align(merge_blocks=True)` or `merge_alignment_blocks()` |
#' | Evaluation | `alignment_accuracy()`, `alignment_f1()`, etc. |
#' | N-dim alignment | Automatic for N<=4 via `align()` |
#' | Progressive | Automatic for N>4 via `align()` |
#' | Prior-guided learning | `bootstrap_matrix(prior_matrix=...)` |
#' | Block-aware learning | `bootstrap_matrix(block_merge=True)` |
