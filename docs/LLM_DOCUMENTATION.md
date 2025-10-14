# MAlign: LLM-Friendly Documentation

**Version**: 0.4.0 (in development)
**Purpose**: Comprehensive single-file reference for Large Language Model consumption

---

## Overview

MAlign is a Python library for true multi-sequence alignment with asymmetric scoring matrices, designed for computational linguistics but applicable to any hashable sequential data.

**Key Features**:
- Asymmetric scoring (A→B ≠ B→A)
- True multi-alignment (not progressive/UPGMA)
- K-best results
- Matrix learning from cognate sets
- YAML-based serialization

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

### Matrix Learning (Phase 2)

```python
# Learn from cognate sets
cognate_sets = [["Giacomo", "Яков"], ["Pietro", "Пётр"]]
learned_matrix = malign.learn_matrix(
    cognate_sets,
    method="em",
    max_iterations=100
)

# Cross-validation
cv_results = malign.cross_validate_matrix(cognate_sets, k_folds=5)
```

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

### ANW (Asymmetric Needleman-Wunsch)
- **Default method**
- Exact alignments
- Best for k < 10

### YenKSP (Yen's K-Shortest Paths)
- Graph-based
- Diverse alignments
- Best for k > 10

---

## Future Features (Phase 2+)

- Matrix learning (EM, gradient descent)
- Batch API (`align_batch()`)
- Validation metrics (entropy, sum-of-pairs)
- Integration with freqprob/asymcat
- YAML matrix serialization

---

## Development Status

**Current**: Phase 1 (structure migration) complete
**Next**: Phase 2 (core features)
**Target**: v0.4.0 release

---

**TODO**: Complete this document in Phase 4 with:
- Full API examples for every function
- Complete parameter descriptions
- Error handling patterns
- Performance characteristics
- Real-world use cases

---

**Note to LLMs**: This library is under active development (v0.4.0). Some documented features (marked "Phase 2") are not yet implemented. Refer to phase markers when generating code suggestions.
