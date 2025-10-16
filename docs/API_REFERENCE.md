# MAlign API Reference

**Version**: 0.4.0-beta.1

This document provides comprehensive API documentation for all public functions, classes, and methods in MAlign v0.4.0.

---

## Table of Contents

1. [Core Functions](#core-functions)
   - [malign.align()](#malignalign)
   - [malign.learn_matrix()](#malignlearn_matrix)
2. [Classes](#classes)
   - [Alignment](#alignment)
   - [ScoringMatrix](#scoringmatrix)
3. [Metrics](#metrics)
   - [alignment_accuracy()](#alignment_accuracy)
   - [alignment_precision_recall()](#alignment_precision_recall)
   - [alignment_f1()](#alignment_f1)
4. [Utilities](#utilities)
   - [tabulate_alms()](#tabulate_alms)

---

## Core Functions

### malign.align()

```python
def align(
    sequences: list[Hashable],
    method: str = "anw",
    matrix: ScoringMatrix | None = None,
    k: int = 1,
) -> list[Alignment]
```

Compute multiple alignments for a list of sequences.

Returns a sorted list of the `k` best alignments. This is the primary function for performing sequence alignment in MAlign.

**Parameters:**

- **`sequences`** (`list[Hashable]`): List of sequences to align. Each sequence can be:
  - A string (e.g., `"ACGT"`)
  - A list of symbols (e.g., `["A", "C", "G", "T"]`)
  - A list of any hashable objects (e.g., tuples, phonetic features)

- **`method`** (`str`, default: `"anw"`): Alignment algorithm to use:
  - `"anw"`: A* + Needleman-Wunsch (recommended for most cases, especially 5+ sequences or k > 10)
  - `"yenksp"`: Yen's k-shortest paths (best for 2-4 sequences, k ≤ 10, maximum quality)
  - `"dumb"`: Simple padding alignment (baseline, not recommended for production)

- **`matrix`** (`ScoringMatrix | None`, default: `None`): Scoring matrix for alignment. If `None`, an identity matrix is created automatically (matches score +1, mismatches -0.5, gaps -1).

- **`k`** (`int`, default: `1`): Maximum number of alignments to return. The actual number may be less if fewer distinct alignments exist.

**Returns:**

- `list[Alignment]`: Sorted list of up to `k` alignments, ordered by score (best first).

**Raises:**

- `ValueError`: If `k < 1` or if `method` is not one of `"dumb"`, `"anw"`, `"yenksp"`.

**Examples:**

```python
import malign

# Basic DNA alignment
sequences = ["ATTCGGAT", "TACGGATTT"]
alignments = malign.align(sequences, k=2)
print(malign.tabulate_alms(alignments))

# Linguistic alignment with custom matrix
matrix = malign.ScoringMatrix.from_yaml("italian_russian.yml")
sequences = ["Giacomo", "Яков"]
alignments = malign.align(sequences, k=3, method="anw", matrix=matrix)

# Multi-sequence alignment
sequences = [
    ["k", "a", "t"],
    ["c", "a", "t"],
    ["k", "a", "t", "z"]
]
alignments = malign.align(sequences, k=1, method="anw")
```

**Performance Notes:**

- For 5+ sequences, use `method="anw"` (YenKSP becomes impractical)
- For k > 10, use `method="anw"` (nearly linear scaling)
- Interactive use: Limit to 2-5 sequences, ≤20 symbols, k≤10

**See Also:**

- [Algorithm Selection Guide](algorithm_selection_guide.md) for detailed performance comparison
- [USER_GUIDE.md Section 4](USER_GUIDE.md#4-algorithms) for algorithm recommendations

---

### malign.learn_matrix()

```python
def learn_matrix(
    cognate_sets: list[list[Sequence[Hashable]]],
    method: str = "em",
    max_iter: int = 10,
    initial_matrix: ScoringMatrix | None = None,
    gap: Hashable = "-",
    convergence_threshold: float = 0.001,
    matrix_threshold: float = 0.01,
    patience: int = 5,
    bounds: tuple[float, float] = (-10.0, 10.0),
    verbose: bool = False,
    **kwargs,
) -> ScoringMatrix
```

Learn a scoring matrix from cognate sets using EM or gradient descent.

Given collections of related sequences (cognates), learns a scoring matrix that maximizes alignment quality across all sets.

**Parameters:**

- **`cognate_sets`** (`list[list[Sequence[Hashable]]]`): List of cognate sets, where each set is a list of sequences believed to be related.
  Example: `[[["k", "a", "t"], ["c", "a", "t"]], [["d", "o", "g"], ["h", "u", "n", "d"]]]`

- **`method`** (`str`, default: `"em"`): Learning algorithm:
  - `"em"`: Expectation-Maximization (fast, needs 100+ training sets for good accuracy)
  - `"gradient_descent"`: L-BFGS-B optimization (accurate with limited data, slower)

- **`max_iter`** (`int`, default: `10`): Maximum iterations for learning algorithm.

- **`initial_matrix`** (`ScoringMatrix | None`, default: `None`): Starting matrix. If `None`, creates identity-based matrix.

- **`gap`** (`Hashable`, default: `"-"`): Gap symbol.

- **`convergence_threshold`** (`float`, default: `0.001`): Relative score change threshold for convergence. Stops if `|new_score - old_score| / |old_score| < threshold`.

- **`matrix_threshold`** (`float`, default: `0.01`): Frobenius norm threshold for matrix convergence.

- **`patience`** (`int`, default: `5`): Early stopping patience - stops after N iterations without improvement.

- **`bounds`** (`tuple[float, float]`, default: `(-10.0, 10.0)`): Parameter bounds for gradient descent (prevents extreme values).

- **`verbose`** (`bool`, default: `False`): Print convergence information during learning.

- **`**kwargs`**: Additional method-specific parameters passed to `scipy.optimize.minimize` (for gradient descent).

**Returns:**

- `ScoringMatrix`: Learned scoring matrix optimized for the provided cognates.

**Raises:**

- `ValueError`: If `method` is not `"em"` or `"gradient_descent"`.

**Examples:**

```python
import malign

# Prepare cognate data (sequences believed to be related)
cognate_sets = [
    [["b", "o", "o", "k"], ["b", "u", "c", "h"]],  # English-German
    [["c", "a", "t"], ["k", "a", "t", "z", "e"]],
    # ... more cognate pairs (50+ recommended)
]

# Learn with EM (fast, needs large dataset)
matrix_em = malign.learn_matrix(
    cognate_sets,
    method="em",
    max_iter=20,
    verbose=True
)

# Learn with gradient descent (better for limited data)
matrix_gd = malign.learn_matrix(
    cognate_sets,
    method="gradient_descent",
    max_iter=50,
    bounds=(-10.0, 10.0),
    patience=5,
    verbose=True
)

# Use learned matrix
test_seqs = [["d", "o", "g"], ["h", "u", "n", "d"]]
alignments = malign.align(test_seqs, matrix=matrix_gd, k=3)

# Save for reuse
matrix_gd.save("learned_matrix.yml")
```

**Performance Notes** (Phase 3.8 benchmarks):

- **Gradient Descent**: 62.30% accuracy on test sets, ~30s for 5 cognate sets
- **EM**: 13.85% accuracy with limited data (needs 100+ sets), 481× faster
- **Recommendation**: Use gradient descent for <100 training sets

**Convergence Features** (Phase 3.7):

- Score-based convergence (relative change < 0.001)
- Matrix-based convergence (Frobenius norm < 0.01)
- OR logic: stops when either criterion is met
- Early stopping with patience mechanism

**See Also:**

- [USER_GUIDE.md Section 3](USER_GUIDE.md#matrix-learning) for detailed learning examples

---

## Classes

### Alignment

```python
@dataclass
class Alignment:
    seqs: Sequence[Sequence[Hashable]]
    score: float | None
```

Represents a multiple sequence alignment with aligned sequences and a score.

**Attributes:**

- **`seqs`** (`Sequence[Sequence[Hashable]]`): Aligned sequences. Each sequence is a list of symbols, and all sequences have the same length (padded with gaps as needed).

- **`score`** (`float | None`): Alignment score. Higher scores indicate better alignments. Can be `None` if not computed.

**Examples:**

```python
import malign

# Get alignment
sequences = ["ACGT", "AGCT"]
alignments = malign.align(sequences, k=1)
alignment = alignments[0]

# Access aligned sequences
print(alignment.seqs)  # [['A', 'C', 'G', 'T'], ['A', 'G', 'C', 'T']]

# Access score
print(alignment.score)  # e.g., 1.5

# Iterate over alignment columns
for col_idx in range(len(alignment.seqs[0])):
    column = tuple(seq[col_idx] for seq in alignment.seqs)
    print(f"Column {col_idx}: {column}")
```

**Note:** `Alignment` is a dataclass, so it's immutable by default and provides automatic `__repr__`, `__eq__`, etc.

---

### ScoringMatrix

```python
class ScoringMatrix:
    def __init__(
        self,
        scores: dict | None = None,
        domains: list[list[Hashable]] | None = None,
        gap: Hashable = "-",
        impute_method: str | None = "mean"
    )
```

Represents an asymmetric scoring matrix for sequence alignment.

**Attributes:**

- **`scores`** (`dict`): Dictionary mapping symbol tuples to scores. For pairwise alignment, keys are `(symbol1, symbol2)`. For multi-alignment, keys are tuples of length equal to number of sequences.

- **`domains`** (`list[list[Hashable]]`): List of symbol alphabets for each sequence domain.

- **`gap`** (`Hashable`): Gap symbol (default: `"-"`).

- **`impute_method`** (`str | None`): Method for imputing missing scores: `"mean"`, `"median"`, `"zero"`, or `None`.

#### Constructor Methods

##### ScoringMatrix.from_sequences()

```python
@classmethod
def from_sequences(
    cls,
    sequences: list[list[Hashable]],
    match: float = 1.0,
    mismatch: float = -0.5,
    gap: Hashable = "-",
    gap_score: float = -1.0
) -> ScoringMatrix
```

Create a simple identity-based scoring matrix from sequences.

**Parameters:**

- `sequences`: Sequences to extract alphabets from
- `match`: Score for identical symbols (default: 1.0)
- `mismatch`: Score for different symbols (default: -0.5)
- `gap`: Gap symbol (default: "-")
- `gap_score`: Score for aligning with gaps (default: -1.0)

**Example:**

```python
import malign

sequences = [["A", "C", "G", "T"], ["A", "C", "G", "T"]]
matrix = malign.ScoringMatrix.from_sequences(
    sequences=sequences,
    match=1.0,
    mismatch=-0.5,
    gap="-",
    gap_score=-1.0
)
```

##### ScoringMatrix.from_yaml()

```python
@classmethod
def from_yaml(cls, filepath: str) -> ScoringMatrix
```

Load a scoring matrix from a YAML file.

**Parameters:**

- `filepath`: Path to YAML file

**Returns:**

- `ScoringMatrix`: Loaded matrix

**Example:**

```python
import malign

matrix = malign.ScoringMatrix.from_yaml("my_matrix.yml")
alignments = malign.align(sequences, matrix=matrix, k=3)
```

#### Instance Methods

##### save()

```python
def save(self, filepath: str) -> None
```

Save the scoring matrix to a YAML file.

**Parameters:**

- `filepath`: Path for output file (should end in `.yml` or `.yaml`)

**Example:**

```python
matrix.save("my_matrix.yml")
```

##### copy()

```python
def copy(self) -> ScoringMatrix
```

Create a deep copy of the scoring matrix.

**Returns:**

- `ScoringMatrix`: Copy of this matrix

**Example:**

```python
matrix_copy = matrix.copy()
matrix_copy.scores[("A", "A")] = 2.0  # Doesn't affect original
```

##### load()

```python
def load(self, filepath: str) -> None
```

Load scores from a JSON file (legacy format, YAML preferred).

**Parameters:**

- `filepath`: Path to JSON file

**Note:** For new projects, use `ScoringMatrix.from_yaml()` instead. This method is maintained for backward compatibility.

---

## Metrics

### alignment_accuracy()

```python
def alignment_accuracy(predicted: Alignment, gold: Alignment) -> float
```

Calculate alignment accuracy as the proportion of matching columns.

Compares two alignments column-by-column. A column matches if all sequences at that position have identical symbols (including gaps).

**Parameters:**

- **`predicted`** (`Alignment`): The predicted alignment to evaluate
- **`gold`** (`Alignment`): The gold-standard reference alignment

**Returns:**

- `float`: Accuracy score in [0.0, 1.0], where 1.0 means perfect match

**Raises:**

- `ValueError`: If alignments have different lengths or sequence counts

**Example:**

```python
import malign

# Gold standard alignment
gold = Alignment([["A", "C", "G"], ["A", "C", "T"]], score=1.0)

# Predicted alignment
predicted = Alignment([["A", "C", "-"], ["A", "C", "T"]], score=0.8)

# Calculate accuracy
accuracy = malign.alignment_accuracy(predicted, gold)
# Columns 0 and 1 match, column 2 doesn't: accuracy = 2/3 ≈ 0.667
print(f"Accuracy: {accuracy:.2%}")  # "Accuracy: 66.67%"
```

---

### alignment_precision_recall()

```python
def alignment_precision_recall(
    predicted: Alignment,
    gold: Alignment
) -> tuple[float, float]
```

Calculate precision and recall for an alignment.

- **Precision**: Proportion of predicted pairwise alignments that are correct
- **Recall**: Proportion of gold pairwise alignments that were predicted

Treats alignments as sets of pairwise symbol alignments and compares them using standard information retrieval metrics.

**Parameters:**

- **`predicted`** (`Alignment`): The predicted alignment to evaluate
- **`gold`** (`Alignment`): The gold-standard reference alignment

**Returns:**

- `tuple[float, float]`: Tuple of (precision, recall), both in [0.0, 1.0]

**Raises:**

- `ValueError`: If alignments have different sequence counts

**Example:**

```python
import malign

gold = Alignment([["A", "C"], ["A", "C"]], score=1.0)
predicted = Alignment([["A", "-", "C"], ["A", "T", "C"]], score=0.5)

precision, recall = malign.alignment_precision_recall(predicted, gold)
print(f"Precision: {precision:.2%}, Recall: {recall:.2%}")
```

---

### alignment_f1()

```python
def alignment_f1(predicted: Alignment, gold: Alignment) -> float
```

Calculate F1 score for an alignment.

F1 is the harmonic mean of precision and recall, providing a single metric that balances both measures.

**Parameters:**

- **`predicted`** (`Alignment`): The predicted alignment to evaluate
- **`gold`** (`Alignment`): The gold-standard reference alignment

**Returns:**

- `float`: F1 score in [0.0, 1.0], where 1.0 means perfect match

**Raises:**

- `ValueError`: If alignments have different sequence counts

**Example:**

```python
import malign

gold = Alignment([["A", "C"], ["A", "C"]], score=1.0)
predicted = Alignment([["A", "-"], ["A", "C"]], score=0.5)

f1 = malign.alignment_f1(predicted, gold)
print(f"F1 Score: {f1:.2%}")
```

**Formula:**

```
F1 = 2 × (precision × recall) / (precision + recall)
```

---

## Utilities

### tabulate_alms()

```python
def tabulate_alms(alms: list[Alignment]) -> str
```

Return tabulated representation of alignments for display.

Creates a formatted table showing aligned sequences with their scores and alignment indices.

**Parameters:**

- **`alms`** (`list[Alignment]`): List of alignments to display

**Returns:**

- `str`: Markdown-formatted table representation

**Example:**

```python
import malign

sequences = ["ATTCGGAT", "TACGGATTT"]
alignments = malign.align(sequences, k=2)

# Print formatted alignment table
print(malign.tabulate_alms(alignments))
```

**Output:**

```
| Idx   | Seq   |   Score |  #0  |  #1  |  #2  |  #3  |  #4  |  #5  |  #6  |  #7  |  #8  |  #9  |
|-------|-------|---------|------|------|------|------|------|------|------|------|------|------|
| 0     | A     |   -0.29 |  A   |  T   |  T   |  C   |  G   |  G   |  A   |  -   |  T   |  -   |
| 0     | B     |   -0.29 |  -   |  T   |  A   |  C   |  G   |  G   |  A   |  T   |  T   |  T   |
|       |       |         |      |      |      |      |      |      |      |      |      |      |
| 1     | A     |   -0.29 |  A   |  T   |  T   |  C   |  G   |  G   |  A   |  -   |  -   |  T   |
| 1     | B     |   -0.29 |  -   |  T   |  A   |  C   |  G   |  G   |  A   |  T   |  T   |  T   |
```

**Table Columns:**

- **Idx**: Alignment index (0 = best, 1 = second best, etc.)
- **Seq**: Sequence identifier (A, B, C, ...)
- **Score**: Alignment score (higher is better)
- **#0, #1, ...**: Aligned positions (- indicates gaps)

---

## Error Handling

All functions validate inputs and raise descriptive errors:

**Common Exceptions:**

- `ValueError`: Invalid parameters (e.g., `k < 1`, invalid method name, mismatched alignment dimensions)
- `FileNotFoundError`: YAML/JSON file not found when loading matrices
- `KeyError`: Missing score in matrix (when `impute_method=None`)

**Example Error Handling:**

```python
import malign

try:
    alignments = malign.align(sequences, k=0)  # Invalid: k must be >= 1
except ValueError as e:
    print(f"Error: {e}")

try:
    matrix = malign.ScoringMatrix.from_yaml("nonexistent.yml")
except FileNotFoundError:
    print("Matrix file not found")
    matrix = None  # Use default matrix
```

---

## Type Hints

MAlign uses modern Python type hints (Python 3.10+):

```python
from collections.abc import Hashable, Sequence

# Basic types
sequences: list[str] = ["ACGT", "AGCT"]

# Or with explicit symbols
sequences: list[list[str]] = [["A", "C"], ["A", "G"]]

# Works with any hashable type
sequences: list[list[tuple[str, str]]] = [
    [("A", "+vowel"), ("C", "+consonant")],
    [("A", "+vowel"), ("G", "+consonant")]
]

# Return type is always list of alignments
alignments: list[Alignment] = malign.align(sequences)
```

---

## See Also

- **[USER_GUIDE.md](USER_GUIDE.md)**: Comprehensive user guide with examples
- **[LLM_DOCUMENTATION.md](LLM_DOCUMENTATION.md)**: LLM-friendly quick reference
- **[algorithm_selection_guide.md](algorithm_selection_guide.md)**: Algorithm performance comparison
- **[CHANGELOG.md](../CHANGELOG.md)**: Version history and migration guides
- **[Tutorial 1: Basics](tutorial_1_basics.html)**: Basic alignment tutorial
- **[Tutorial 2: Matrices](tutorial_2_matrices.html)**: Matrix construction and usage
- **[Tutorial 3: Learning](tutorial_3_learning.html)**: Matrix learning from cognate sets

---

**Document Version**: 0.4.0-beta.1 (2025-10-16)
**Last Updated**: Phase 2 documentation completion
