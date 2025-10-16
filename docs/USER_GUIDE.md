# MAlign User Guide

**Version**: 0.4.0-beta.1

This guide covers the essential features and usage patterns for MAlign v0.4.0. For detailed API documentation, see [API_REFERENCE.md](API_REFERENCE.md).

---

## Table of Contents

1. [Getting Started](#getting-started)
2. [Core Concepts](#core-concepts)
3. [Matrix Management](#matrix-management)
4. [Algorithms](#algorithms)
5. [Advanced Topics](#advanced-topics)
6. [Performance & Limits](#performance--limits)
7. [Troubleshooting](#troubleshooting)
8. [API Reference](#api-reference)

---

## 1. Getting Started

### Installation

```bash
pip install malign
```

For development or beta versions:

```bash
pip install malign==0.4.0b1
```

### Quick Start

The simplest use case is aligning two or more sequences:

```python
import malign

# Basic DNA sequence alignment
sequences = ["ATTCG", "ATGC"]
alignments = malign.align(sequences, k=2)

# View results in tabular format
print(malign.tabulate_alms(alignments))
```

Output:
```
| Idx   | Seq   |   Score |  #0  |  #1  |  #2  |  #3  |  #4  |  #5  |
|-------|-------|---------|------|------|------|------|------|------|
| 0     | A     |   -0.45 |  A   |  T   |  T   |  C   |  G   |  -   |
| 0     | B     |   -0.45 |  A   |  T   |  -   |  -   |  G   |  C   |
```

### Linguistic Example

MAlign is particularly useful for linguistic applications, comparing words across languages:

```python
import malign

# Compare Italian and Russian cognates
sequences = ["Giacomo", "Яков"]  # Both mean "Jacob"

# Align with multiple alternatives
alignments = malign.align(sequences, k=3, method="anw")

# View all alternatives
for i, aln in enumerate(alignments):
    print(f"\nAlignment {i+1} (score: {aln.score:.2f}):")
    print(malign.tabulate_alms([aln]))
```

### Multiple Sequences

Align three or more sequences simultaneously (not progressive alignment):

```python
import malign

# Align three cognates at once
sequences = [
    ["k", "a", "t"],      # English "cat"
    ["c", "a", "t"],      # Latin "cattus"
    ["k", "a", "t", "z"]  # German "Katze"
]

# Get best alignment
alignments = malign.align(sequences, k=1, method="anw")
print(malign.tabulate_alms(alignments))
```

### Key Parameters

- **`sequences`**: List of sequences (strings or lists of tokens)
- **`k`**: Number of best alignments to return (default: 1)
- **`method`**: Algorithm to use - `"anw"` (default) or `"yenksp"`
- **`matrix`**: Scoring matrix (optional, defaults to identity-based scoring)

### Working with Custom Alphabets

MAlign works with any hashable Python objects:

```python
import malign

# Phonetic features
seq1 = [("p", "+stop", "+bilabial"), ("a", "+vowel", "+open")]
seq2 = [("b", "+stop", "+bilabial"), ("a", "+vowel", "+open")]

alignments = malign.align([seq1, seq2], k=1)
```

---

## 2. Core Concepts

### Asymmetric Scoring

Unlike traditional alignment algorithms (e.g., BLAST, ClustalW), MAlign supports **asymmetric scoring matrices** where the score for aligning symbol A to symbol B can differ from aligning B to A.

**Why Asymmetry Matters:**

In linguistic sound change, transformations are often directional:
- Latin *p* → Spanish *b* (voicing) is common
- Spanish *b* → Latin *p* (devoicing) is rare

A symmetric matrix cannot capture this directionality. MAlign allows you to specify:
```
score(p, b) = 0.8   # Common change
score(b, p) = -0.2  # Rare/unlikely change
```

**Domain-Specific Alignment:**

MAlign can align sequences from different domains (e.g., IPA transcriptions vs. orthography):

```python
import malign

# Define asymmetric matrix
matrix = malign.ScoringMatrix()
matrix.scores = {
    ("th", "θ"): 2.0,   # English orthography → IPA
    ("θ", "th"): 0.5,   # IPA → orthography (less natural)
    ("t", "t"): 1.0,    # Identity match
    # ... more mappings
}

# Align orthographic and phonetic forms
seqs = [["th", "i", "n", "k"], ["θ", "ɪ", "ŋ", "k"]]
alignments = malign.align(seqs, matrix=matrix, k=1)
```

### True Multi-Alignment

MAlign performs **true multiple sequence alignment**, not progressive pairwise alignment.

**Comparison with Progressive Methods:**

| Feature | MAlign | Progressive (UPGMA, ClustalW) |
|---------|--------|-------------------------------|
| **Method** | Simultaneous multi-way alignment | Series of pairwise alignments |
| **Scoring** | Global probability across all sequences | Sum of pairwise scores |
| **Optimality** | Globally optimal (with ANW/YenKSP) | Dependent on guide tree |
| **Asymmetry** | Full support | Not supported |
| **Use Case** | Linguistic cognates, directional changes | Biological sequences |

**Why It Matters:**

Progressive alignment can introduce artifacts where early decisions constrain later alignments. MAlign evaluates all sequences together, finding alignments that maximize the overall probability/score.

```python
import malign

# Three related sequences
sequences = [
    ["k", "a", "t"],
    ["c", "a", "t"],
    ["k", "a", "t", "z"]
]

# MAlign finds best global alignment
# (not "align seq1-seq2, then add seq3")
alignments = malign.align(sequences, k=1, method="anw")
```

### K-Best Results

MAlign can return multiple high-scoring alignments (k-best solutions), useful for:

**1. Exploring Alternatives:**
When the optimal alignment is ambiguous, examine multiple possibilities:

```python
# Get top 5 alignments
alignments = malign.align(sequences, k=5)

for i, aln in enumerate(alignments):
    print(f"Alternative {i+1} (score: {aln.score:.2f})")
    print(malign.tabulate_alms([aln]))
```

**2. Uncertainty Quantification:**
If the top k alignments have similar scores, the alignment is uncertain:

```python
alignments = malign.align(sequences, k=3)

scores = [aln.score for aln in alignments]
score_range = max(scores) - min(scores)

if score_range < 0.1:
    print("High uncertainty: multiple alignments are nearly equally good")
```

**3. Downstream Analysis:**
Use multiple alignments for consensus building or phylogenetic reconstruction:

```python
# Extract all possible alignment columns
all_columns = []
for aln in malign.align(sequences, k=10):
    for col_idx in range(len(aln.seqs[0])):
        column = tuple(seq[col_idx] for seq in aln.seqs)
        all_columns.append(column)

# Analyze column frequencies
from collections import Counter
freq = Counter(all_columns)
```

---

## 3. Matrix Management

### Construction Methods

#### 1. From Sequences (Simple Identity Matrix)

The fastest way to create a matrix for quick testing:

```python
import malign

# Create identity-based matrix from sequences
sequences = [["A", "C", "G", "T"], ["A", "C", "G", "T"]]

matrix = malign.ScoringMatrix.from_sequences(
    sequences=sequences,
    match=1.0,      # Score for identical symbols
    mismatch=-0.5,  # Score for different symbols
    gap="-",
    gap_score=-1.0  # Score for gap alignment
)

# Use in alignment
alignments = malign.align(sequences, matrix=matrix, k=1)
```

#### 2. Manual Construction

For complete control over all alignment scores:

```python
import malign

# Define domains (symbol alphabets for each sequence)
domains = [
    ["-", "A", "C", "G", "T"],  # Sequence 1 symbols
    ["-", "A", "C", "G", "T"]   # Sequence 2 symbols
]

# Create empty matrix
matrix = malign.ScoringMatrix(domains=domains, gap="-")

# Manually set scores
matrix.scores[("A", "A")] = 2.0   # Match
matrix.scores[("A", "C")] = -1.0  # Mismatch
matrix.scores[("A", "-")] = -2.0  # Gap penalty
# ... set all other scores
```

#### 3. From Pre-computed Probabilities

If you have probability data, convert to log-odds scores:

```python
import numpy as np
import malign

# Your probability matrix (e.g., from frequency counts)
prob_matrix = {
    ("A", "A"): 0.25,
    ("A", "C"): 0.05,
    # ...
}

# Convert to log-odds scores
matrix = malign.ScoringMatrix(domains=domains, gap="-")
for (sym1, sym2), prob in prob_matrix.items():
    # Log-odds: log(P(align) / P(random))
    background_prob = 0.1  # Adjust based on your data
    matrix.scores[(sym1, sym2)] = np.log(prob / background_prob)
```

#### 4. From Cognate Sets (Learning)

Learn optimal matrices from training data (see [Matrix Learning](#matrix-learning) below):

```python
import malign

# Training data: cognate sets
cognate_sets = [
    [["k", "a", "t"], ["c", "a", "t"]],
    [["d", "o", "g"], ["d", "o", "c", "k"]],
    # ... more cognate pairs
]

# Learn matrix from data
matrix = malign.learn_matrix(
    cognate_sets,
    method="gradient_descent",
    max_iter=50
)
```

### Serialization (YAML Format)

**Saving Matrices:**

```python
import malign

# Create or learn a matrix
matrix = malign.ScoringMatrix.from_sequences(
    sequences=[["A", "C", "G"], ["A", "C", "G"]],
    match=1.0,
    mismatch=-0.5,
    gap="-",
    gap_score=-1.0
)

# Save to YAML (primary format)
matrix.save("my_matrix.yml")
```

**Loading Matrices:**

```python
import malign

# Load from YAML
matrix = malign.ScoringMatrix.from_yaml("my_matrix.yml")

# Use immediately
alignments = malign.align(sequences, matrix=matrix, k=3)
```

**YAML Format Structure:**

```yaml
domains:
  - ["-", "A", "C", "G", "T"]
  - ["-", "A", "C", "G", "T"]
gap: "-"
scores:
  "('-', '-')": -1.0
  "('A', 'A')": 1.0
  "('A', 'C')": -0.5
  "('A', '-')": -1.0
  # ... all score pairs
```

**Backward Compatibility:**

JSON format is still supported for backward compatibility:

```python
# Load old JSON matrices
matrix = malign.ScoringMatrix()
matrix.load("old_matrix.json")

# Convert to YAML
matrix.save("new_matrix.yml")
```

### Matrix Imputation

Handle missing scores with automatic imputation:

```python
import malign

# Create matrix with incomplete scores
matrix = malign.ScoringMatrix(
    domains=[["-", "A", "C"], ["-", "A", "C"]],
    gap="-",
    impute_method="mean"  # or "median", "zero"
)

# Only set some scores
matrix.scores[("A", "A")] = 1.0
matrix.scores[("C", "C")] = 1.0

# Missing scores automatically imputed when used
# (impute_method controls how missing values are filled)
```

### Matrix Learning

**New in v0.4.0**: Learn scoring matrices from cognate data.

#### Expectation-Maximization (EM)

Fast iterative learning, good for initial matrices:

```python
import malign

# Cognate sets: sequences believed to be related
cognate_sets = [
    [["b", "o", "o", "k"], ["b", "u", "c", "h"]],  # English-German
    [["c", "a", "t"], ["k", "a", "t", "z", "e"]],   # English-German
    # ... 50+ more pairs recommended
]

# Learn with EM
matrix = malign.learn_matrix(
    cognate_sets,
    method="em",
    max_iter=20,
    verbose=True  # Show convergence info
)

# Use learned matrix
test_seqs = [["d", "o", "g"], ["h", "u", "n", "d"]]
alignments = malign.align(test_seqs, matrix=matrix, k=3)
```

#### Gradient Descent

More accurate but slower, finds global optimum:

```python
import malign

# Learn with gradient descent (L-BFGS-B)
matrix = malign.learn_matrix(
    cognate_sets,
    method="gradient_descent",
    max_iter=50,
    bounds=(-10.0, 10.0),  # Prevent extreme values
    patience=5,            # Early stopping
    verbose=True
)
```

**Learning Parameters:**

- `method`: `"em"` (fast) or `"gradient_descent"` (accurate)
- `max_iter`: Maximum iterations (default: 10)
- `convergence_threshold`: Score change threshold (default: 0.001)
- `patience`: Early stopping after N iterations without improvement (default: 5)
- `verbose`: Print convergence information (default: False)

**Performance Note** (from Phase 3.8 benchmarks):
- **Gradient Descent**: 62.30% accuracy on test sets, ~30s for 5 cognate sets
- **EM**: 13.85% accuracy with limited data (needs 100+ training sets), 481× faster

**Recommendation**: Use gradient descent for better accuracy when training data is limited (<100 sets). Use EM for large datasets (100+ sets) where speed matters.

---

## 4. Algorithms

MAlign provides two alignment algorithms with different trade-offs. For detailed guidance, see [Algorithm Selection Guide](algorithm_selection_guide.md).

### ANW (A* + Needleman-Wunsch)

**Best for:** General-purpose alignment, especially with 5+ sequences or large k values.

```python
import malign

sequences = [["A", "C", "G"], ["A", "G", "C"], ["A", "C", "T"]]
alignments = malign.align(sequences, k=10, method="anw")
```

**Characteristics:**
- **Speed**: Faster, especially for k > 10
- **Quality**: High quality, near-optimal alignments
- **Scaling**: Nearly linear with k (~2× slower for k=20 vs k=1)
- **Use Case**: Large problems, interactive applications, production systems

**Performance** (from Phase 3 benchmarks):
- 2 sequences: 0.001s
- 5 sequences: 0.4s (~400× slower)
- k=1 to k=20: ~2× slower (nearly linear)
- Practical limits: 5-8 sequences, k < 50

### YenKSP (Yen's k-Shortest Paths)

**Best for:** Maximum quality on small problems (2-4 sequences, k ≤ 10).

```python
import malign

sequences = [["k", "a", "t"], ["c", "a", "t"]]
alignments = malign.align(sequences, k=5, method="yenksp")
```

**Characteristics:**
- **Speed**: Slower, especially for many sequences
- **Quality**: Highest quality, exhaustive search with theoretical guarantees
- **Scaling**: Exponential with sequence count (362× for 5 vs 2 sequences)
- **Use Case**: Research, benchmarking, gold standards, batch processing

**Performance** (from Phase 3 benchmarks):
- 2 sequences: Very fast
- 5 sequences: 362× slower than 2 sequences
- Practical limits: ≤ 4 sequences, k ≤ 10

### Quick Decision Guide

```
How many sequences?
├─ 2-4 sequences
│  ├─ k ≤ 10? → Use YenKSP (maximum quality)
│  └─ k > 10? → Use ANW (faster for large k)
└─ 5+ sequences → Use ANW (YenKSP too slow)
```

**Examples:**

```python
import malign

# Small problem, want maximum quality
sequences = [["cat"], ["bat"]]
alms = malign.align(sequences, k=5, method="yenksp")

# Large problem or many alternatives
sequences = [["s1"], ["s2"], ["s3"], ["s4"], ["s5"]]
alms = malign.align(sequences, k=20, method="anw")

# Interactive use (need fast response)
alms = malign.align(sequences, k=3, method="anw")
```

**For complete algorithm comparison including:**
- Detailed performance benchmarks
- Sequence length and k value scaling
- Production vs research recommendations
- Complete decision flowchart

**See:** [Algorithm Selection Guide](algorithm_selection_guide.md)

---

## 5. Advanced Topics

### Batch Processing

**TODO**: Document align_batch() in Phase 4.

### Validation Metrics

**TODO**: Document all metrics in Phase 4.

### Integration with External Libraries

- **freqprob**: TODO
- **asymcat**: TODO
- **nhandu**: TODO

---

## 6. Performance & Limits

### Computational Limits

Based on Phase 3 benchmarking results:

**Recommended Configuration (ANW method):**
- **Sequence count**: 2-5 sequences (practical), 6-8 sequences (batch only)
- **Sequence length**: Up to 20 symbols (interactive), 50 symbols (batch)
- **k value**: k ≤ 10 (interactive), k ≤ 50 (batch)

**Maximum Practical Limits:**
- **Sequences**: 8 sequences possible but slow (minutes per alignment)
- **Length**: 100+ symbols possible but memory-intensive
- **k value**: k = 100+ possible with ANW, limited by memory

### Scaling Characteristics

**Sequence Count** (ANW, k=1):
```
2 sequences:  0.001s  (baseline)
3 sequences:  0.004s  (4× slower)
4 sequences:  0.053s  (53× slower)
5 sequences:  0.397s  (397× slower)
```
→ **Exponential scaling**: Each additional sequence multiplies time by ~7-10×

**Sequence Length** (ANW, 3 sequences, k=1):
```
5 symbols:   0.004s  (baseline)
10 symbols:  0.004s  (same)
15 symbols:  0.005s  (1.25× slower)
20 symbols:  0.010s  (2.5× slower)
```
→ **Sub-quadratic scaling**: Doubling length ~2.5× slower

**K Value** (ANW, 3 sequences):
```
k=1:   0.004s  (baseline)
k=5:   0.004s  (same)
k=10:  0.006s  (1.5× slower)
k=20:  0.008s  (2× slower)
```
→ **Nearly linear scaling**: Practical up to k=50

### Performance Tips

**1. Choose the Right Algorithm:**
```python
# For many sequences or large k, use ANW
if num_sequences >= 5 or k > 10:
    method = "anw"
else:
    method = "yenksp"  # Higher quality for small problems
```

**2. Limit k Value:**
```python
# Don't compute more alignments than needed
alignments = malign.align(sequences, k=5)  # Not k=100
```

**3. Batch Processing:**
```python
# Process multiple alignment tasks in parallel
from multiprocessing import Pool

def align_task(sequences):
    return malign.align(sequences, k=3, method="anw")

with Pool(4) as pool:
    results = pool.map(align_task, sequence_batches)
```

**4. Optimize Matrices:**
```python
# Learned matrices often produce better (faster) alignments
# than default identity matrices
matrix = malign.learn_matrix(training_data, method="gradient_descent")
alignments = malign.align(sequences, matrix=matrix, k=3)
```

**5. Profile Your Use Case:**
```python
import time

start = time.time()
alignments = malign.align(sequences, k=10, method="anw")
elapsed = time.time() - start

print(f"Aligned {len(sequences)} sequences in {elapsed:.3f}s")
```

### Benchmarks Summary

From `scripts/benchmarks.py` (Phase 3):

| Configuration | ANW Time | YenKSP Time | Recommendation |
|---------------|----------|-------------|----------------|
| 2 seqs, k=1   | 0.001s   | 0.001s      | Either method OK |
| 3 seqs, k=1   | 0.004s   | 0.008s      | Either method OK |
| 4 seqs, k=1   | 0.053s   | 0.12s       | Either method OK |
| 5 seqs, k=1   | 0.397s   | 7.22s       | **Use ANW** |
| 3 seqs, k=10  | 0.006s   | 0.015s      | Either method OK |
| 3 seqs, k=20  | 0.008s   | 0.035s      | **Use ANW for k>10** |

**Key Insight**: YenKSP provides marginally better quality but becomes impractical for 5+ sequences or k > 10. ANW is the recommended default for most use cases.

---

## 7. Troubleshooting

### Common Issues

**TODO**: Populate based on user feedback.

### Performance Issues

**TODO**: Add profiling and optimization tips.

---

## 8. API Reference

See [API_REFERENCE.md](API_REFERENCE.md) for complete API documentation.

---

**Document Status**: Essential sections completed for v0.4.0-beta.1. Advanced topics (sections 5, 7) to be expanded based on user feedback.
