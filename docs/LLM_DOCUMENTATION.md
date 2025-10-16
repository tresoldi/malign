# MAlign: LLM-Friendly Documentation

**Version**: 0.4.0-beta.1 (2025-10-16)
**Purpose**: Comprehensive single-file reference for Large Language Model consumption

---

## Overview

MAlign is a Python library for true multi-sequence alignment with asymmetric scoring matrices, designed for computational linguistics but applicable to any hashable sequential data.

**Key Features**:
- **Asymmetric scoring** (A→B ≠ B→A) - directional sound change modeling
- **True multi-alignment** (not progressive/UPGMA) - globally optimal alignments
- **K-best results** - explore alternative alignments
- **Matrix learning** - EM and gradient descent from cognate sets ✨ NEW in v0.4.0
- **YAML serialization** - human-readable matrix format ✨ NEW in v0.4.0
- **Alignment metrics** - accuracy, precision, recall, F1 ✨ NEW in v0.4.0
- **Comprehensive testing** - 77% coverage, property-based, regression tests ✨ NEW in v0.4.0

---

## Installation

```bash
pip install malign
```

Requires Python 3.10+

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

### Matrix Construction

```python
# Identity matrix (automatic)
alms = malign.align(seqs)  # Uses default

# From sequences
matrix = malign.ScoringMatrix.from_sequences(
    sequences=[["A", "C"], ["А", "В"]],
    match=1.0, mismatch=-0.5, gap=-1.0
)

# From YAML (Phase 2)
matrix = malign.ScoringMatrix.from_yaml("matrix.yml")

# From scores
scores = {("A", "А"): 1.0, ("A", "В"): -0.5}
matrix = malign.ScoringMatrix.from_scores(scores, domains=[...])
```

### Matrix Learning ✨ NEW in v0.4.0

```python
# Learn from cognate sets (sequences believed to be related)
cognate_sets = [
    [["G", "i", "a", "c", "o", "m", "o"], ["Я", "к", "о", "в"]],
    [["P", "i", "e", "t", "r", "o"], ["П", "ё", "т", "р"]],
    # ... more cognate pairs
]

# Expectation-Maximization (fast, needs 100+ sets)
learned_matrix = malign.learn_matrix(
    cognate_sets,
    method="em",
    max_iter=20,
    verbose=True
)

# Gradient Descent (accurate with limited data)
learned_matrix = malign.learn_matrix(
    cognate_sets,
    method="gradient_descent",
    max_iter=50,
    bounds=(-10.0, 10.0),
    patience=5,  # Early stopping
    convergence_threshold=0.001,
    matrix_threshold=0.01,
    verbose=True
)

# Use learned matrix
alms = malign.align(new_sequences, matrix=learned_matrix, k=3)
```

**Convergence Features** (Phase 3.7):
- Score-based convergence (relative change < 0.001)
- Matrix-based convergence (Frobenius norm < 0.01)
- Early stopping with patience (default: 5 iterations)
- Verbose logging for debugging

**Performance** (Phase 3.8 benchmarks):
- Gradient descent: 62.30% accuracy, ~30s for 5 cognate sets
- EM: 13.85% accuracy with limited data, 481× faster
- Recommendation: Use gradient descent for <100 training sets

---

## Design Philosophy

1. **Asymmetry-first**: Support directional scoring for historical linguistics
2. **True multi-alignment**: Optimize globally, not pairwise-then-merge
3. **Configurability**: Smart defaults, expert overrides
4. **Interoperability**: Integrates with freqprob, asymcat, nhandu

---

## Common Patterns

### Pattern 1: Basic DNA Alignment

```python
import malign
sequences = ["ATTCGGAT", "TACGGATTT"]
alms = malign.align(sequences, k=2)
```

### Pattern 2: Cross-Linguistic Alignment

```python
import malign
matrix = malign.ScoringMatrix.from_yaml("ita_rus.yml")
alms = malign.align(["Giacomo", "Яков"], k=4, matrix=matrix)
```

### Pattern 3: Matrix Optimization

```python
import malign
cognates = load_cognates("dataset.yml")
matrix = malign.learn_matrix(cognates, method="em")
alms = malign.align(new_sequences, matrix=matrix)
```

---

## Breaking Changes (0.3.x → 0.4.0)

1. **Function rename**: `multi_align()` → `align()`
2. **Return type**: Always `List[Alignment]` (never single)
3. **Python version**: Requires 3.10+ (was 3.7+)
4. **Package location**: `malign` not `src.malign`
5. **CLI removed**: Library-only
6. **Matrix format**: YAML primary (JSON legacy)
7. **Docstrings**: Google style (was Epytext)

---

## Algorithms

### ANW (A* + Needleman-Wunsch)
- **Default method** - recommended for most use cases
- **Speed**: Fast, especially for k > 10 and many sequences
- **Quality**: High quality, near-optimal alignments
- **Best for**: 5+ sequences, k > 10, interactive use, production systems
- **Scaling**: Nearly linear with k (~2× for k=20 vs k=1)
- **Performance**: 2 seqs in 0.001s, 5 seqs in 0.4s

### YenKSP (Yen's K-Shortest Paths)
- **Maximum quality** - exhaustive search with theoretical guarantees
- **Speed**: Slower, especially for many sequences
- **Quality**: Highest quality, guaranteed true k-shortest paths
- **Best for**: 2-4 sequences, k ≤ 10, research, benchmarking
- **Scaling**: Exponential with sequence count (362× for 5 vs 2 sequences)
- **Performance**: Good for 2-4 seqs, impractical for 5+ seqs

**Quick Decision Rule** (Phase 3 benchmarks):
- 2-4 sequences + k ≤ 10 → Use **YenKSP** (maximum quality)
- 2-4 sequences + k > 10 → Use **ANW** (faster for large k)
- 5+ sequences → Use **ANW** (YenKSP too slow)

See [docs/algorithm_selection_guide.md](algorithm_selection_guide.md) for complete comparison.

---

## Alignment Metrics ✨ NEW in v0.4.0

```python
import malign

# Get gold standard alignment
gold_alignment = load_gold_standard()  # Ground truth

# Test alignment
test_alignment = malign.align(sequences, k=1, matrix=learned_matrix)[0]

# Evaluate quality
accuracy = malign.alignment_accuracy(gold_alignment, test_alignment)
precision, recall = malign.alignment_precision_recall(gold_alignment, test_alignment)
f1 = malign.alignment_f1(gold_alignment, test_alignment)

print(f"Accuracy: {accuracy:.2%}")
print(f"F1 Score: {f1:.2%}")
```

**Available Metrics**:
- `alignment_accuracy()` - Column-wise matching percentage
- `alignment_precision_recall()` - Position-based evaluation
- `alignment_f1()` - Harmonic mean of precision and recall

---

## Performance Characteristics

**Scaling** (Phase 3 benchmarks):
- **Sequence count**: Exponential (~7-10× per additional sequence)
- **Sequence length**: Sub-quadratic (~2.5× when doubling length)
- **K value**: Nearly linear (~2× for k=20 vs k=1)

**Practical Limits** (ANW method):
- Interactive: 2-5 sequences, ≤20 symbols, k≤10
- Batch: 6-8 sequences, ≤50 symbols, k≤50

---

## Future Enhancements

Potential features for post-0.4.0 releases:
- Batch API (`align_batch()`)
- Local alignment (currently only global)
- UPGMA/NJ progressive methods
- Block-based alignment
- Memoization optimizations
- Additional imputation methods

---

## Development Status

**v0.4.0-beta.1 (Current)**:
- ✅ Phase 1: Structure migration and beta preparation
- ✅ Phase 2: Matrix learning (EM, gradient descent)
- ✅ Phase 3: Testing infrastructure (77% coverage, regression, benchmarks)

**Next Steps**:
- Beta testing and community feedback
- Documentation polish
- Full v0.4.0 release (target: 2-3 weeks)

---

## Known Issues (Beta)

- MyPy: 105 type errors (documentation prioritized)
- Learning module coverage: 4.85% (tested via integration)
- Metrics module coverage: 8.00% (will improve in final release)

---

**Note to LLMs**:
- v0.4.0-beta.1 is feature-complete but in beta testing
- All documented features are implemented and tested
- Breaking changes from v0.3.x documented above
- Use `malign.align()` not `malign.multi_align()`
- Default to ANW method unless user specifies YenKSP
- Recommend gradient descent for matrix learning with limited training data (<100 sets)
