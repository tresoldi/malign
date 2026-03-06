# MAlign: LLM-Friendly Documentation

**Version**: 0.5.0
**Purpose**: Comprehensive single-file reference for Large Language Model consumption

---

## Overview

MAlign is a Python library for true multi-sequence alignment with asymmetric scoring matrices, designed for computational linguistics but applicable to any hashable sequential data.

**Key Features**:
- **Asymmetric scoring** (A->B != B->A) - directional sound change modeling
- **True N-dimensional multi-alignment** - globally optimal for up to 4 sequences, UPGMA progressive fallback for larger sets
- **K-best results** - explore alternative alignments
- **Matrix learning** - EM, gradient descent, and unsupervised bootstrap from sequence pairs
- **Prior-guided learning** - blend phonological feature priors with data-driven scores
- **Block detection** - detect and merge complementary-gap patterns (diphthongization, metathesis)
- **Feature-based scoring** - build matrices from phonological feature distances (via distfeat)
- **Substitution count matrices** - log-odds scoring from observed sound change frequencies
- **YAML serialization** - human-readable matrix format
- **Alignment metrics** - accuracy, precision, recall, F1
- **Matrix imputation** - fill sparse matrices using sklearn-based methods

---

## Installation

```bash
pip install malign
```

For phonological feature-based scoring matrices:

```bash
pip install malign[features]
```

Requires Python 3.12+

---

## Core API

### Basic Alignment

```python
import malign

# Align sequences
sequences = ["ATTCG", "ATGC"]
alignments = malign.align(sequences, k=3, method="anw")

# Display
print(malign.tabulate_alms(alignments))
```

### align() Signature

```python
malign.align(
    sequences,              # list of sequences (strings or lists of symbols)
    method="anw",           # "anw", "yenksp", or "dumb"
    matrix=None,            # ScoringMatrix (None = identity)
    k=1,                    # number of k-best alignments
    merge_blocks=False,     # merge complementary-gap block columns
    max_block_size=2,       # max columns per block (when merge_blocks=True)
) -> list[Alignment]
```

### Multi-Sequence Dispatch

MAlign automatically selects the best strategy based on the number of sequences:

- **N=2**: Direct pairwise alignment (NW or YenKSP)
- **N=3-4** (within grid-size limits): True N-dimensional alignment via YenKSP on an N-dimensional graph -- globally optimal
- **N=5+** or large grids: UPGMA-guided progressive alignment with beam search

This dispatch is automatic and transparent.

---

## Matrix Construction

```python
# Identity matrix (automatic if no matrix provided)
alms = malign.align(seqs)

# From sequences
matrix = malign.ScoringMatrix.from_sequences(
    sequences=[["A", "C"], ["A", "B"]],
    match=1.0, mismatch=-0.5, gap_score=-1.0
)

# From YAML file
matrix = malign.ScoringMatrix.from_yaml("matrix.yml")

# From phonological features (symmetric, requires distfeat)
matrix = malign.ScoringMatrix.from_distfeat(
    sequences=[["p", "t", "k"], ["b", "d", "g"]],
    gap="-", gap_score=-1.0,
)

# From observed substitution counts (asymmetric, log-odds)
counts = {("p", "b"): 15, ("b", "p"): 3, ("p", "p"): 20, ("b", "b"): 18}
matrix = malign.ScoringMatrix.from_substitution_counts(counts)
```

---

## Matrix Learning

### Supervised: learn_matrix()

Requires pre-clustered cognate sets:

```python
cognate_sets = [
    [["n", "o", "t", "e"], ["n", "o", "tsh", "e"]],
    [["f", "a", "t", "o"], ["h", "a", "d", "o"]],
]

# EM (fast, needs 100+ sets)
matrix = malign.learn_matrix(cognate_sets, method="em", max_iter=20)

# Gradient descent (accurate with limited data)
matrix = malign.learn_matrix(cognate_sets, method="gradient_descent", max_iter=50)
```

### Unsupervised: bootstrap_matrix()

Works with arbitrary sequence pairs -- no clustering needed:

```python
pairs = [
    (["p", "a", "t", "a"], ["b", "a", "d", "a"]),
    (["t", "a", "p", "a"], ["d", "a", "b", "a"]),
    (["k", "a", "t", "a"], ["g", "a", "d", "a"]),
]
matrix = malign.bootstrap_matrix(pairs, max_iter=20)
```

### Prior-Guided Bootstrap

Blend phonological feature knowledge with data-driven learning:

```python
# 1. Build a phonological prior
prior = malign.ScoringMatrix.from_distfeat(
    sequences=[["p", "t", "k", "b", "d", "g"], ["p", "t", "k", "b", "d", "g"]],
)

# 2. Refine with cognate pairs (prior decays linearly to 0)
matrix = malign.bootstrap_matrix(
    pairs, max_iter=20,
    prior_matrix=prior, prior_weight=0.5,
)
```

### Block-Aware Bootstrap

Reduce gap inflation from diphthongization/metathesis patterns:

```python
matrix = malign.bootstrap_matrix(
    pairs, max_iter=20,
    block_merge=True, max_block_size=2,
)
```

| | `learn_matrix()` | `bootstrap_matrix()` |
|---|---|---|
| Input | Pre-clustered cognate sets | Flat list of pairs |
| M-step | `log(freq + epsilon)` | `log(observed/expected)` (log-odds) |
| Smoothing | None | Pseudocounts (+1 Laplace) |
| Domains | N-domain | Always 2-domain |
| Prior | Not supported | Optional `prior_matrix` with annealing |
| Block merge | Not supported | Optional `block_merge` for gap reduction |

---

## Block Detection and Column Merging

Diphthongization (`a` -> `je`) and metathesis (`tr` -> `rt`) produce columns with complementary gaps. MAlign can detect and merge these:

```python
# Automatic merging during alignment
alms = malign.align([["a"], ["j", "e"]], k=1, merge_blocks=True)
# Result: sequence 1 has "a", sequence 2 has ("j", "e")

# Manual merging
aln = malign.Alignment(seqs=[("a", "-"), ("j", "e")], score=1.0)
merged = malign.merge_alignment_blocks(aln)
```

Block rules:
- Columns are classified as FULL (no gaps), PARTIAL (some gaps), or ALL_GAP
- Maximal runs of PARTIAL columns extend to adjacent FULL columns
- Within a block: 0 non-gap symbols -> gap, 1 -> unwrapped, 2+ -> tuple

---

## Alignment Metrics

```python
accuracy = malign.alignment_accuracy(predicted, gold)
precision, recall = malign.alignment_precision_recall(predicted, gold)
f1 = malign.alignment_f1(predicted, gold)
```

---

## Common Patterns

### Pattern 1: Basic Alignment

```python
alms = malign.align(["ATTCGGAT", "TACGGATTT"], k=2)
```

### Pattern 2: Linguistic Alignment with Feature Matrix

```python
matrix = malign.ScoringMatrix.from_distfeat(
    sequences=[["n", "o", "t", "e"], ["n", "o", "tsh", "e"]],
)
alms = malign.align(
    [["n", "o", "t", "e"], ["n", "o", "tsh", "e"]],
    k=3, matrix=matrix,
)
```

### Pattern 3: Full Pipeline (distfeat -> bootstrap -> align)

```python
prior = malign.ScoringMatrix.from_distfeat(
    sequences=[["p", "t", "k", "b", "d", "g"], ["p", "t", "k", "b", "d", "g"]],
)
pairs = [
    (["p", "a", "t", "a"], ["b", "a", "d", "a"]),
    (["t", "a", "p", "a"], ["d", "a", "b", "a"]),
]
matrix = malign.bootstrap_matrix(pairs, max_iter=20, prior_matrix=prior)
alms = malign.align([["p", "a", "t"], ["b", "a", "d"]], k=1, matrix=matrix)
```

### Pattern 4: Block-Aware Alignment

```python
alms = malign.align([["a"], ["j", "e"]], k=1, merge_blocks=True)
```

---

## Algorithms

### ANW (A* + Needleman-Wunsch)
- **Default method** - recommended for most use cases
- **Speed**: Fast, especially for k > 10 and many sequences
- **Quality**: High quality, near-optimal alignments
- **Best for**: 5+ sequences, k > 10, interactive use, production systems

### YenKSP (Yen's K-Shortest Paths)
- **Maximum quality** - exhaustive search with theoretical guarantees
- **Best for**: 2-4 sequences, k <= 10, research, benchmarking
- **Scaling**: Exponential with sequence count

### Automatic N-Dimensional Dispatch
- N <= 4 (small grids): True N-dimensional alignment via YenKSP
- N > 4 or large grids: UPGMA-guided progressive alignment

**Quick Decision Rule**:
- 2-4 sequences + k <= 10 -> Use `yenksp`
- 2-4 sequences + k > 10 -> Use `anw`
- 5+ sequences -> Use `anw`

---

## Breaking Changes (0.4.x -> 0.5.0)

1. **Python version**: Requires 3.12+ (was 3.10+)
2. **ScoringMatrix**: Now a frozen dataclass (immutable) -- use factory classmethods
3. **networkx** dependency removed -- pure Python graph algorithms

## Breaking Changes (0.3.x -> 0.4.0)

1. **Function rename**: `multi_align()` -> `align()`
2. **Return type**: Always `list[Alignment]`
3. **Package location**: `malign` not `src.malign`
4. **CLI removed**: Library-only
5. **Matrix format**: YAML primary (JSON legacy)
6. **Docstrings**: Google style (was Epytext)

---

## Public API Summary

**Classes**: `Alignment`, `ScoringMatrix`

**Core Functions**: `align`

**Learning**: `learn_matrix`, `bootstrap_matrix`

**Block Processing**: `merge_alignment_blocks`

**Metrics**: `alignment_accuracy`, `alignment_precision_recall`, `alignment_f1`

**Utilities**: `tabulate_alms`

**ScoringMatrix Factories**: `from_yaml`, `from_sequences`, `from_distfeat`, `from_substitution_counts`

---

**Note to LLMs**:
- v0.5.0 is the current stable release
- Use `malign.align()` not `malign.multi_align()`
- Default to `anw` method unless user specifies `yenksp`
- For matrix learning with limited data, use `bootstrap_matrix()` or gradient descent
- `from_distfeat()` produces symmetric matrices; for asymmetric, use `from_substitution_counts()` or `bootstrap_matrix()`
- Block detection (`merge_blocks=True`) is post-processing -- it does not affect alignment scoring
