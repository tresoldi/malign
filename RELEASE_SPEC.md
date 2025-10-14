# MAlign 0.4.0 Release Specification

**Date**: 2025-10-14
**Target Version**: 0.4.0
**Timeline**: Flexible (side project, 3-4 months focused work)
**Breaking Changes**: Acceptable and expected

---

## Executive Summary

MAlign is a Python library for true multi-sequence alignment with asymmetric scoring matrices, designed primarily for computational linguistics but applicable to any sequential hashable data. Version 0.4.0 represents a major refactoring and enhancement for public release and eventual JOSS submission.

**Key Changes in 0.4.0**:
- Modern Python (3.10+) with strict type checking
- Flat package structure following PROJECT_STRUCTURE_TEMPLATE.md
- YAML-based matrix serialization (human-readable, version-control friendly)
- Matrix learning/optimization from cognate sets
- Google-style docstrings
- Builder pattern for API construction
- Comprehensive testing (>90% coverage)
- Nhandu-based tutorial documentation
- Library-only (CLI removed)

---

## Core Value Propositions

1. **Asymmetric Scoring**: Support for directional measures (A→B ≠ B→A), essential for:
   - Historical sound change modeling (directionality in evolution)
   - Cross-alphabet alignment (e.g., Latin ↔ Cyrillic have inherent asymmetry)
   - Any domain where transitions have inherent directionality

2. **True Multi-Sequence Alignment**: Unlike progressive methods (UPGMA, NJ), MAlign:
   - Considers all pairwise alignments simultaneously
   - Optimizes global alignment score across all sequences
   - Returns k-best alternatives for inspection
   - Provides better accuracy and interpretability for linguistic data

3. **Matrix Learning from Data**: Given cognate sets or parallel sequences:
   - Learn optimal scoring matrices via EM or gradient descent
   - Joint optimization of alignments and scores
   - Maximize alignment scores across assumed-correct cognate sets
   - Integration with freqprob and asymcat libraries

---

## Target Users & Use Cases

### Primary Users
- **Computational linguists** doing historical/comparative linguistics
- **Researchers** needing asymmetric sequence alignment (any domain with hashable sequences)

### Typical Workloads
- **Small**: Single collection of ~20 sequences (most common)
- **Medium**: ~100 collections × ~100 sequences each (intensive linguistics work)
- **Not supported**: Massive genomic-scale datasets (out of scope)

### Usage Patterns
- Manual inspection of k-best alignments
- Semi-automated analysis of cognate sets
- Matrix optimization from gold-standard data
- Experimentation with different scoring approaches

---

## Resolved Decisions

### Project Structure
- **✓ Package layout**: Flat structure at root (`malign/`, not `src/malign/`)
- **✓ CLI**: Complete removal (library-only)
- **✓ Python version**: 3.10+ (breaking change from 3.7+)
- **✓ Docstrings**: Migrate to Google style (from Epytext)
- **✓ Build system**: pyproject.toml (replace setup.py)

### API Design
- **✓ Main function**: Rename `multi_align` → `align` (clearer, shorter)
- **✓ Return type**: Always `List[Alignment]` even for k=1 (consistency)
- **✓ Matrix construction**: Builder pattern (`ScoringMatrix.from_*()` methods)
- **✓ Batch API**: Add `align_batch()` for multiple alignment jobs
- **✓ Validation metrics**: Both as Alignment methods and standalone functions

### Matrix Serialization
- **✓ Primary format**: YAML (human-readable, version-control friendly)
- **✓ Legacy support**: Keep JSON loading for backward compatibility
- **✓ Serialization**: Support both YAML and JSON export

### Code Quality
- **✓ Linting**: ruff (format + check)
- **✓ Type checking**: mypy with strict settings
- **✓ Coverage**: >90% required (enforced in CI)
- **✓ Testing**: pytest with markers (slow, integration, etc.)

### Documentation
- **✓ Platform**: Nhandu-based tutorials (no ReadTheDocs, no mkdocs)
- **✓ Format**: `docs/tutorial_*.py` → HTML via nhandu
- **✓ Examples**: Use nhandu package for realistic data

### Refactoring Approach
- **✓ ScoringMatrix**: Keep monolithic for 0.4.0, extract `matrix_learning.py`
- **✓ Validation metrics**: Both Alignment methods + utils functions
- **✓ UPGMA/NJ**: Defer to 0.5.0+ (focus on core differentiators)

---

## Technical Architecture

### Algorithm Portfolio

#### 1. Asymmetric Needleman-Wunsch (`anw`, default)
- **Best for**: Exact alignments, smaller k values
- **Characteristics**: More conservative, proven correctness
- **Status**: Mature, needs performance benchmarking

#### 2. Yen's K-Shortest Paths (`yenksp`)
- **Best for**: Larger k values, exploring alternatives
- **Characteristics**: Graph-based, more diverse paths
- **Status**: Functional, needs usage guidance

#### 3. Dumb Baseline (`dumb`)
- **Purpose**: Simple baseline for benchmarking
- **Status**: Keep for testing/comparison

#### 4. UPGMA/NJ (deferred to 0.5.0)
- **Purpose**: Fast approximate method for 10+ sequences
- **Positioning**: Complement, not replacement
- **Status**: Not in 0.4.0 scope

### Performance Analysis (Phase 3 Task)
**Action Required**: Empirical analysis to provide users with:
- Complexity analysis (time/space) for ANW vs YenKSP
- Performance benchmarks on linguistic data
- Decision flowchart in documentation
- Computational limits (maximum sequences: 5? 10? 20?)

---

## Project Structure Migration

### Target Structure (Following PROJECT_STRUCTURE_TEMPLATE.md)

```
malign/                          # Main package (flat at root)
├── __init__.py                  # Version + exports
├── alignment.py                 # Alignment dataclass + validation methods
├── anw.py                       # ANW algorithm
├── dumb.py                      # Baseline method
├── malign.py                    # Main entry (align, align_batch)
├── matrix_learning.py           # NEW - Matrix optimization/learning
├── scoring_matrix.py            # Core matrix class (refactored)
├── utils.py                     # Utilities + validation functions
├── yenksp.py                    # YenKSP algorithm
└── py.typed                     # NEW - PEP 561 type marker

docs/
├── tutorial_1_basics.py         # NEW - Nhandu tutorial
├── tutorial_2_matrices.py       # NEW - Matrix construction
├── tutorial_3_learning.py       # NEW - Matrix learning
├── tutorial_4_advanced.py       # NEW - Advanced usage
├── figures/                     # NEW - Tutorial figures
├── API_REFERENCE.md             # NEW - Complete API docs
├── USER_GUIDE.md                # NEW - Comprehensive guide
└── LLM_DOCUMENTATION.md         # NEW - LLM-friendly reference

tests/
├── test_alignment.py            # NEW - Alignment class tests
├── test_anw.py                  # NEW - ANW algorithm
├── test_dumb.py                 # NEW - Baseline
├── test_integration.py          # NEW - Full pipelines
├── test_malign.py               # KEEP - Main functions (update)
├── test_matrix.py               # RENAME from test_matrix.py
├── test_matrix_learning.py      # NEW - Optimization tests
├── test_property_based.py       # NEW - Hypothesis tests
├── test_regression.py           # NEW - Linguistic gold data
├── test_utils.py                # NEW - Utility functions
├── test_yenksp.py               # NEW - YenKSP algorithm
├── data/                        # NEW - Test fixtures/data
└── htmlcov/                     # Coverage reports (gitignored)

scripts/
└── benchmarks.py                # NEW - Performance benchmarks

.github/
├── workflows/
│   └── quality.yml              # UPDATE - CI/CD pipeline
└── dependabot.yml               # NEW - Dependency updates

# Root files
├── Makefile                     # NEW - Development commands
├── pyproject.toml               # NEW - Replace setup.py + requirements.txt
├── CHANGELOG.md                 # NEW - Keep a Changelog format
├── README.md                    # UPDATE - Comprehensive rewrite
├── CONTRIBUTING.md              # KEEP - Update for new structure
├── LICENSE                      # KEEP - MIT
├── .gitignore                   # UPDATE - Comprehensive Python ignore
└── RELEASE_SPEC.md              # KEEP - This document

# Files to DELETE
├── old_demo.py                  # DELETE
├── workbench.py                 # DELETE
├── setup.py                     # DELETE (→ pyproject.toml)
├── requirements.txt             # DELETE (→ pyproject.toml)
├── MANIFEST.in                  # DELETE (not needed)
├── src/                         # DELETE (move to flat malign/)
└── src/malign/__main__.py       # DELETE (no CLI)
```

### Migration Checklist (Part of Phase 1)

- [ ] Create new `malign/` directory at root
- [ ] Move all modules from `src/malign/` to `malign/`
- [ ] Delete `src/` directory
- [ ] Delete `__main__.py` (CLI removal)
- [ ] Create `malign/py.typed` (empty file for PEP 561)
- [ ] Create `Makefile` following template
- [ ] Create `pyproject.toml` (replace setup.py)
- [ ] Migrate dependencies from requirements.txt to pyproject.toml
- [ ] Delete setup.py, requirements.txt, MANIFEST.in
- [ ] Update .gitignore comprehensively
- [ ] Create docs/ structure with tutorial templates
- [ ] Create scripts/benchmarks.py stub
- [ ] Update GitHub workflows
- [ ] Create CHANGELOG.md
- [ ] Delete old_demo.py and workbench.py

---

## Matrix Serialization Format

### YAML Format Specification (Primary)

**Design Goals**:
- Human-readable and editable
- Version-control friendly (line-based diffs)
- Supports multi-domain matrices
- Backward compatible (can load old JSON)

#### Format Structure

```yaml
# MAlign Scoring Matrix
# Format version: 0.4.0
version: "0.4.0"
gap: "-"
impute_method: mean  # or null for no imputation
num_domains: 2

domains:
  - [A, C, G, T, "-"]
  - [А, В, Г, Т, "-"]

scores:
  # Format: [symbol_0, symbol_1, ..., score]
  # All symbols must be from their respective domains
  - [A, А, 1.0]
  - [A, В, -0.5]
  - [A, Г, -0.3]
  - [A, Т, -0.4]
  - [A, "-", -1.0]
  # ... more scores

  # Can use YAML anchors for repeated values
  - [C, В, -0.5]
  - [C, Г, 0.8]
  # ... etc

# Optional metadata
metadata:
  description: "Italian-Russian phonetic alignment matrix"
  created_by: "malign.optimize_matrix"
  created_date: "2025-10-14"
  training_data: "cognate_sets_ita_rus.yml"
```

#### Sub-Matrix Representation

For pairwise sub-matrices (domains not fully specified):

```yaml
version: "0.4.0"
gap: "-"
num_domains: 2

domains:
  - [A, C, G, T, "-"]
  - [А, В, Г, Т, "-"]

scores:
  # Standard pairwise scores
  - [A, А, 1.0]
  - [A, В, -0.5]
  # ...
```

#### Multi-Domain (3+) Matrices

```yaml
version: "0.4.0"
gap: "-"
num_domains: 3

domains:
  - [a, b, c, "-"]  # Sequence 1 alphabet
  - [x, y, z, "-"]  # Sequence 2 alphabet
  - [α, β, γ, "-"]  # Sequence 3 alphabet

scores:
  # Three-way alignment scores
  - [a, x, α, 2.5]
  - [a, x, β, 1.0]
  # Can include partial (sub-matrix) specifications with null
  - [a, x, null, 0.5]  # Score for (a,x,?) any third symbol
  # ...
```

### JSON Format (Legacy Support)

Keep existing JSON format for backward compatibility:

```json
{
  "version": "0.4.0",
  "gap": "-",
  "domain_range": [0, 1],
  "domains": [
    ["A", "C", "G", "T", "-"],
    ["А", "В", "Г", "Т", "-"]
  ],
  "scores": {
    "A / А": 1.0,
    "A / В": -0.5,
    "A / NULL": -1.0
  }
}
```

### API for Serialization

```python
# Loading
matrix = ScoringMatrix.from_file("matrix.yml")  # Auto-detect format
matrix = ScoringMatrix.from_yaml("matrix.yml")
matrix = ScoringMatrix.from_json("matrix.json")  # Legacy

# Saving
matrix.save("matrix.yml")  # Default to YAML
matrix.save("matrix.yml", format="yaml")
matrix.save("matrix.json", format="json")  # Legacy export
```

---

## API Design (v0.4.0)

### Core Alignment API

```python
import malign

# Basic usage (RENAMED from multi_align)
alignments = malign.align(
    ["Giacomo", "Яков"],
    k=4,
    method="anw",  # or "yenksp", "dumb"
    matrix=scorer
)
# Always returns List[Alignment], even if k=1

# Tabulation (unchanged)
print(malign.tabulate_alms(alignments))

# Batch processing (NEW)
results = malign.align_batch(
    sequences=[
        ["Giacomo", "Яков"],
        ["Pietro", "Пётр"],
        ["Giovanni", "Иван"]
    ],
    k=5,
    method="anw",
    matrix=scorer
)
# Returns List[List[Alignment]]
```

### Alignment Object (Enhanced)

```python
alm = alignments[0]

# Existing attributes
alm.seqs       # Tuple[Sequence[Hashable], ...]
alm.score      # float
len(alm)       # Number of sequences
alm[0]         # Get sequence by index

# NEW: Validation metrics (both methods and standalone functions)
alm.sum_of_pairs_score(matrix)      # float
alm.entropy()                         # List[float] - per site
alm.consensus(threshold=0.5)          # Sequence[Hashable]
alm.to_dict()                         # Dict representation

# Alternatively, use standalone functions
score = malign.sum_of_pairs(alm, matrix)
entropy = malign.alignment_entropy(alm)
consensus = malign.alignment_consensus(alm, threshold=0.5)
```

### Matrix Construction (Builder Pattern)

```python
from malign import ScoringMatrix

# From file (auto-detect format)
matrix = ScoringMatrix.from_file("matrix.yml")

# From YAML (explicit)
matrix = ScoringMatrix.from_yaml("matrix.yml")

# From JSON (legacy)
matrix = ScoringMatrix.from_json("matrix.json")

# From sequences (identity matrix)
matrix = ScoringMatrix.from_sequences(
    sequences=[["A", "T", "G"], ["А", "Т", "Г"]],
    match=1.0,
    mismatch=-0.5,
    gap=-1.0
)

# From scores (programmatic construction)
scores = {
    ("A", "А"): 1.0,
    ("A", "В"): -0.5,
    # ...
}
matrix = ScoringMatrix.from_scores(
    scores=scores,
    domains=[["A", "C", "G", "T"], ["А", "В", "Г", "Т"]],
    gap="-",
    impute_method="mean"  # or None
)

# From freqprob distribution (NEW)
from freqprob import FreqProbDist
fp_dist = FreqProbDist.from_sequences(all_sequences)
matrix = ScoringMatrix.from_freqprob(
    dist=fp_dist,
    sequences=sequences,
    gap="-"
)

# From asymcat analysis (NEW)
from asymcat import AsymmetricCategorical
ac = AsymmetricCategorical.from_data(transitions)
matrix = ScoringMatrix.from_asymcat(
    model=ac,
    domains=domains,
    gap="-"
)

# Configurable unobserved symbol handling (NEW)
matrix = ScoringMatrix.from_scores(
    scores=scores,
    domains=domains,
    unobserved_strategy="min"  # "min", "impute", "error"
)
```

### Matrix Learning (NEW Core Feature)

```python
import malign

# Basic matrix optimization from cognate sets
optimized_matrix = malign.learn_matrix(
    cognate_sets=[
        ["Giacomo", "Яков"],
        ["Pietro", "Пётр"],
        ["Giovanni", "Иван"],
        # ... more cognate sets
    ],
    method="em",                  # or "gradient_descent"
    initial_matrix=None,          # Optional starting point
    max_iterations=100,
    convergence_threshold=0.001,
    gap="-"
)

# With external priors
from freqprob import FreqProbDist
fp_dist = FreqProbDist.from_sequences(all_sequences)

matrix = malign.learn_matrix(
    cognate_sets=cognate_sets,
    method="em",
    prior=fp_dist,                # Bayesian prior
    max_iterations=50
)

# Cross-validation for evaluation
cv_results = malign.cross_validate_matrix(
    cognate_sets=cognate_sets,
    k_folds=5,
    method="em",
    max_iterations=100
)
# Returns: {
#   "mean_score": float,
#   "std_score": float,
#   "fold_scores": List[float],
#   "best_matrix": ScoringMatrix
# }
```

### Breaking Changes Summary

1. **Function rename**: `multi_align()` → `align()`
2. **Return type**: Always `List[Alignment]` (never single Alignment)
3. **Matrix construction**: Direct `ScoringMatrix()` discouraged, use builders
4. **Python version**: Requires 3.10+ (dropped 3.7-3.9 support)
5. **Package location**: `malign` not `src.malign`
6. **CLI removed**: No `__main__.py` or command-line interface
7. **Docstring style**: Google (migrated from Epytext)
8. **Matrix format**: YAML primary (JSON still supported)

---

## Matrix Learning & Optimization

### Philosophy
**Core Assumption**: If provided cognate sets are correct, the alignment scores across all sets should be collectively maximized.

### Approach

#### 1. Expectation-Maximization (EM)
**Priority**: Implement first (simpler, more robust)

**Algorithm**:
```
Initialize: Start with identity matrix or user-provided matrix
Repeat until convergence:
  E-step: Align all cognate sets with current matrix
  M-step: Update matrix scores to maximize total alignment score
          (frequency-based: co-aligned symbols get higher scores)
Check: Convergence if score change < threshold
```

**Implementation location**: `malign/matrix_learning.py`

#### 2. Gradient Descent
**Priority**: Implement after EM (more flexible, but complex)

**Algorithm**:
```
Initialize: Random or identity matrix
Define loss: L = -Σ(alignment_scores across all cognate sets)
Repeat until convergence:
  Compute gradient: ∂L/∂matrix_params
  Update: matrix_params -= learning_rate * gradient
  Project: Ensure scores stay in valid range
Check: Convergence if loss change < threshold or gradient norm small
```

**Challenges**:
- Alignment scoring may not be differentiable (use approximation)
- Constrain matrix values to reasonable ranges
- Learning rate tuning

#### 3. Integration with External Tools

**freqprob** (frequency/probability distributions):
```python
from freqprob import FreqProbDist

# Use as Bayesian prior
fp_dist = FreqProbDist.from_sequences(all_sequences)
matrix = malign.learn_matrix(
    cognate_sets,
    method="em",
    prior=fp_dist,  # Regularize with frequency-based prior
)
```

**asymcat** (asymmetric categorical analysis):
```python
from asymcat import AsymmetricCategorical

# Learn directional transition model, use as initialization
ac = AsymmetricCategorical.from_data(observed_transitions)
initial_matrix = ScoringMatrix.from_asymcat(ac, domains, gap="-")

matrix = malign.learn_matrix(
    cognate_sets,
    method="em",
    initial_matrix=initial_matrix,  # Start from asymcat model
)
```

### Validation

```python
# Cross-validation
cv_results = malign.cross_validate_matrix(
    cognate_sets=gold_cognates,
    k_folds=5,
    method="em"
)

# Manual split
train_sets, test_sets = malign.split_cognates(gold_cognates, test_size=0.2)
learned_matrix = malign.learn_matrix(train_sets, method="em")

# Evaluate on test set
test_score = malign.evaluate_matrix(learned_matrix, test_sets)
```

### Test Data
- **Status**: Maintainer will provide gold-standard cognate sets
- **Location**: `tests/data/cognates/`
- **Format**: YAML files with cognate groups
- **Use**: Regression tests, matrix learning validation, examples

---

## Validation Metrics

### Implementation Strategy
**Approach**: Both Alignment methods (convenience) and standalone functions (reusability)

### Metrics to Implement

#### 1. Sum-of-Pairs Score
**Definition**: Sum of all pairwise alignment scores

```python
# Method
score = alignment.sum_of_pairs_score(matrix)

# Function
score = malign.sum_of_pairs(alignment, matrix)

# Returns: float
```

**Use**: Overall alignment quality measure

#### 2. Per-Site Entropy
**Definition**: Shannon entropy at each alignment position

```python
# Method
entropy = alignment.entropy()

# Function
entropy = malign.alignment_entropy(alignment)

# Returns: List[float] - one per alignment position
```

**Use**: Identify conserved vs variable positions

#### 3. Consensus Sequence
**Definition**: Most common symbol at each position (above threshold)

```python
# Method
consensus = alignment.consensus(threshold=0.5)

# Function
consensus = malign.alignment_consensus(alignment, threshold=0.5)

# Returns: Sequence[Hashable] - consensus sequence
```

**Use**: Derive representative sequence from alignment

#### 4. Agreement with Gold Standard (if available)
**Definition**: Alignment accuracy compared to gold reference

```python
# Function only (requires gold reference)
accuracy = malign.alignment_accuracy(alignment, gold_alignment)

# Returns: float in [0, 1]
```

**Use**: Evaluation against known-correct alignments

#### 5. Gap Statistics
**Definition**: Gap counts, gap percentage, gap distribution

```python
# Method
gap_stats = alignment.gap_statistics()

# Returns: Dict with keys:
# {
#   "total_gaps": int,
#   "gap_percentage": float,
#   "gaps_per_sequence": List[int],
#   "gap_runs": List[int],  # Gap run lengths
# }
```

**Use**: Assess alignment quality (too many gaps = poor alignment)

---

## Testing Strategy

### Coverage Goals
- **Target**: >90% line coverage (enforced in CI)
- **Tools**: pytest, pytest-cov, pytest-xdist (parallel)

### Test Organization

#### Test Files Structure
```
tests/
├── test_alignment.py            # Alignment class + metrics
├── test_anw.py                  # ANW algorithm
├── test_dumb.py                 # Dumb baseline
├── test_integration.py          # Full pipelines (multi-step)
├── test_malign.py               # Main align() and align_batch()
├── test_matrix.py               # ScoringMatrix class
├── test_matrix_learning.py      # EM, gradient descent, CV
├── test_property_based.py       # Hypothesis property tests
├── test_regression.py           # Gold-standard linguistic data
├── test_utils.py                # Utility functions
├── test_yenksp.py               # YenKSP algorithm
└── data/                        # Test fixtures
    ├── matrices/                # Test scoring matrices
    ├── cognates/                # Gold cognate sets
    └── sequences/               # Test sequences
```

#### Test Markers (pytest.ini)
```python
@pytest.mark.slow          # Long-running tests (skip in fast mode)
@pytest.mark.integration   # Multi-component tests
@pytest.mark.regression    # Gold-standard validation tests
```

### Test Categories

#### 1. Unit Tests
- Each module independently
- Edge cases: empty sequences, single sequence, all gaps, unobserved symbols
- Matrix operations: imputation, submatrices, serialization

#### 2. Integration Tests
```python
def test_full_alignment_pipeline():
    """Test: sequences → matrix learning → alignment → validation."""
    # Learn matrix from training cognates
    matrix = malign.learn_matrix(train_cognates, method="em")

    # Align test sequences
    alms = malign.align(test_seqs, k=5, matrix=matrix)

    # Validate
    assert len(alms) == 5
    assert alms[0].sum_of_pairs_score(matrix) > threshold
```

#### 3. Regression Tests (Gold Data)
```python
def test_italian_russian_cognates():
    """Validate against gold-standard Italian-Russian alignments."""
    gold_alignments = load_gold_data("ita_rus_cognates.yml")

    for gold in gold_alignments:
        alms = malign.align(gold.sequences, k=1)
        accuracy = malign.alignment_accuracy(alms[0], gold.alignment)
        assert accuracy > 0.85  # 85% agreement threshold
```

#### 4. Property-Based Tests (Hypothesis)
```python
from hypothesis import given, strategies as st

@given(
    seqs=st.lists(
        st.text(alphabet="ACGT", min_size=1, max_size=20),
        min_size=2,
        max_size=5
    )
)
def test_alignment_length_consistency(seqs):
    """All alignments should have same length across sequences."""
    alms = malign.align(seqs, k=10)

    for alm in alms:
        lengths = [len(seq) for seq in alm.seqs]
        assert len(set(lengths)) == 1  # All same length

@given(
    seqs=st.lists(st.text(alphabet="ACGT", min_size=1), min_size=2)
)
def test_alignment_preserves_symbols(seqs):
    """Alignment should preserve original symbols (ignore gaps)."""
    alms = malign.align(seqs, k=1)

    for idx, seq in enumerate(seqs):
        aligned = [s for s in alms[0].seqs[idx] if s != "-"]
        assert aligned == list(seq)
```

#### 5. Performance Tests
```python
@pytest.mark.slow
def test_large_sequence_alignment():
    """Benchmark alignment performance on large sequences."""
    import time

    seqs = generate_sequences(n=10, length=100)

    start = time.time()
    alms = malign.align(seqs, k=5, method="anw")
    duration = time.time() - start

    assert duration < 60.0  # Should complete in under 1 minute
    assert len(alms) == 5
```

### Validation Data Requirements
- **Gold alignments**: 50+ cognate sets with hand-verified alignments
- **Test matrices**: 10+ pre-built matrices for different language pairs
- **Edge cases**: Empty, single-char, all-gap sequences
- **Format**: YAML for readability and version control

---

## Documentation (Nhandu-Based)

### Structure

**No ReadTheDocs, no mkdocs** - Use Nhandu for executable tutorials

```
docs/
├── tutorial_1_basics.py         # Getting started
├── tutorial_2_matrices.py       # Matrix construction and usage
├── tutorial_3_learning.py       # Matrix learning from cognates
├── tutorial_4_advanced.py       # Advanced features (batch, metrics)
├── figures/                     # Auto-generated plots
├── API_REFERENCE.md             # Complete API documentation
├── USER_GUIDE.md                # Comprehensive usage guide
└── LLM_DOCUMENTATION.md         # LLM-friendly reference
```

### Tutorial Content Plan

#### Tutorial 1: Basics (`tutorial_1_basics.py`)
```python
#' # MAlign Tutorial 1: Basic Alignment
#'
#' ## Installation
#'
#' ```bash
#' pip install malign
#' ```

import malign

#' ## Simple Alignment
#'
#' Let's align two DNA sequences:

sequences = ["ATTCG", "ATGC"]
alignments = malign.align(sequences, k=3, method="anw")

#' Visualize results:
print(malign.tabulate_alms(alignments))

#' ## Understanding k-best Alignments
#' ... (more content)
```

#### Tutorial 2: Matrices (`tutorial_2_matrices.py`)
- Matrix construction (identity, from file, from scores)
- YAML format walkthrough
- Imputation methods
- Asymmetric scoring examples

#### Tutorial 3: Learning (`tutorial_3_learning.py`)
- Matrix learning from cognate sets
- EM vs gradient descent
- Cross-validation
- Integration with freqprob/asymcat

#### Tutorial 4: Advanced (`tutorial_4_advanced.py`)
- Batch processing
- Validation metrics
- Performance tuning
- Custom algorithms

### API Reference (docs/API_REFERENCE.md)
Auto-generate from docstrings:
```bash
# Use pydoc-markdown or similar
pydoc-markdown -o docs/API_REFERENCE.md
```

### User Guide (docs/USER_GUIDE.md)
Comprehensive manual covering:
1. **Getting Started**: Installation, quick start
2. **Core Concepts**: Asymmetry, multi-alignment, k-best
3. **Matrix Management**: Construction, serialization, optimization
4. **Algorithms**: ANW vs YenKSP, when to use each
5. **Advanced Topics**: Batch processing, validation, integration
6. **Troubleshooting**: Common issues, performance tips
7. **API Reference**: Link to generated docs

### LLM Documentation (docs/LLM_DOCUMENTATION.md)
Single-file comprehensive reference for LLM consumption:
- Full API specification
- Examples for each function
- Design philosophy
- Common patterns

---

## Code Quality & Standards

### Python Standards (Updated)
- **Python version**: 3.10+ (type hints with modern syntax)
- **Docstrings**: **Google style** (migrated from Epytext)
- **Linting**: ruff (format + check, replaces black/flake8/isort)
- **Type checking**: mypy with strict settings (`disallow_untyped_defs`)
- **Testing**: pytest with markers, >90% coverage enforced
- **Formatting**: ruff format (120 char line length)

### Docstring Migration Example

**Before (Epytext)**:
```python
def align(sequences, k=1, method="anw", matrix=None):
    """
    Compute multiple alignments for sequences.

    @param sequences: List of sequences to align
    @param k: Number of alignments to return
    @param method: Alignment method to use
    @param matrix: Scoring matrix
    @return: List of alignments
    """
```

**After (Google)**:
```python
def align(
    sequences: List[Sequence[Hashable]],
    k: int = 1,
    method: str = "anw",
    matrix: Optional[ScoringMatrix] = None
) -> List[Alignment]:
    """Compute multiple alignments for sequences.

    Args:
        sequences: List of sequences to align (any hashable elements).
        k: Number of best alignments to return. Defaults to 1.
        method: Alignment method ("anw", "yenksp", "dumb"). Defaults to "anw".
        matrix: Scoring matrix. If None, identity matrix is created.

    Returns:
        List of Alignment objects, sorted by score (best first).
        Always returns a list, even if k=1.

    Raises:
        ValueError: If k < 1 or method is invalid.

    Examples:
        >>> alms = malign.align(["ACGT", "AGCT"], k=2)
        >>> print(alms[0].score)
        2.5
    """
```

### Coding Principles (from CLAUDE.md)
1. **Simplicity first (KISS)**
2. **Smart defaults, expert overrides**
3. **Correctness & robustness**
4. **Single responsibility & composability**
5. **Modularity with clean interfaces**
6. **Transparency & clarity**
7. **Representation over logic**
8. **DRY & reuse**
9. **Functional over OO** (prefer pure functions, dataclasses)
10. **Configuration that scales**
11. **Prototype, then optimize**
12. **Interoperability as default**

---

## Phase 1: Detailed Checklist (Structure & Quality)

### 1.1 Project Structure Migration
- [ ] Create `malign/` directory at project root
- [ ] Move all `.py` files from `src/malign/` to `malign/`
- [ ] Delete `src/` directory entirely
- [ ] Delete `malign/__main__.py` (CLI removal)
- [ ] Create empty `malign/py.typed` file (PEP 561 marker)
- [ ] Verify all imports still work after move
- [ ] Update any relative imports if broken

### 1.2 Build System Migration
- [ ] Create `pyproject.toml` from template
- [ ] Migrate metadata from `setup.py` to `pyproject.toml`
- [ ] Migrate dependencies from `requirements.txt` to `pyproject.toml`
- [ ] Configure ruff in `pyproject.toml`
- [ ] Configure mypy in `pyproject.toml`
- [ ] Configure pytest in `pyproject.toml`
- [ ] Configure coverage in `pyproject.toml`
- [ ] Test build: `python -m build`
- [ ] Delete `setup.py`
- [ ] Delete `requirements.txt`
- [ ] Delete `MANIFEST.in`

### 1.3 Development Infrastructure
- [ ] Create `Makefile` from template
- [ ] Customize Makefile for malign-specific needs
- [ ] Create comprehensive `.gitignore` (Python + nhandu)
- [ ] Update `.github/workflows/quality.yml` for new structure
- [ ] Create `.github/dependabot.yml`
- [ ] Test Makefile targets: `make help`, `make quality`, `make test`

### 1.4 Documentation Setup
- [ ] Create `docs/` directory
- [ ] Create `docs/tutorial_1_basics.py` stub
- [ ] Create `docs/tutorial_2_matrices.py` stub
- [ ] Create `docs/tutorial_3_learning.py` stub
- [ ] Create `docs/tutorial_4_advanced.py` stub
- [ ] Create `docs/figures/` directory
- [ ] Create `docs/API_REFERENCE.md` stub
- [ ] Create `docs/USER_GUIDE.md` stub
- [ ] Create `docs/LLM_DOCUMENTATION.md` stub

### 1.5 Test Structure
- [ ] Rename `tests/test_matrix.py` → keep, update
- [ ] Create `tests/test_alignment.py`
- [ ] Create `tests/test_anw.py`
- [ ] Create `tests/test_dumb.py`
- [ ] Create `tests/test_yenksp.py`
- [ ] Create `tests/test_utils.py`
- [ ] Create `tests/test_integration.py`
- [ ] Create `tests/test_matrix_learning.py` (stub)
- [ ] Create `tests/test_property_based.py`
- [ ] Create `tests/test_regression.py` (stub, awaiting data)
- [ ] Create `tests/data/` directory structure
- [ ] Update `tests/test_malign.py` for new API

### 1.6 Code Quality - Ruff
- [ ] Run `ruff check .` and document all issues
- [ ] Fix import sorting issues
- [ ] Fix code style issues (E, W, F rules)
- [ ] Fix complexity issues (C4, SIM rules)
- [ ] Address print statements (T20) - remove or log
- [ ] Fix all other ruff errors
- [ ] Run `ruff format .` to auto-format
- [ ] Verify `make quality` passes ruff checks

### 1.7 Code Quality - MyPy
- [ ] Add `py.typed` marker (done in 1.1)
- [ ] Run `mypy malign/ tests/` and document all errors
- [ ] Add missing type hints to `malign/__init__.py`
- [ ] Add missing type hints to `malign/alignment.py`
- [ ] Add missing type hints to `malign/malign.py`
- [ ] Add missing type hints to `malign/scoring_matrix.py`
- [ ] Add missing type hints to `malign/anw.py`
- [ ] Add missing type hints to `malign/yenksp.py`
- [ ] Add missing type hints to `malign/dumb.py`
- [ ] Add missing type hints to `malign/utils.py`
- [ ] Configure mypy to ignore third-party library issues
- [ ] Verify `make quality` passes mypy checks

### 1.8 Docstring Migration
- [ ] Migrate `malign/__init__.py` to Google style
- [ ] Migrate `malign/alignment.py` to Google style
- [ ] Migrate `malign/malign.py` to Google style
- [ ] Migrate `malign/scoring_matrix.py` to Google style
- [ ] Migrate `malign/anw.py` to Google style
- [ ] Migrate `malign/yenksp.py` to Google style
- [ ] Migrate `malign/dumb.py` to Google style
- [ ] Migrate `malign/utils.py` to Google style
- [ ] Verify all public functions have complete docstrings

### 1.9 Cleanup
- [ ] Delete `old_demo.py`
- [ ] Delete `workbench.py`
- [ ] Remove or update all TODO comments
- [ ] Update `__version__` to "0.4.0.dev0" (development version)
- [ ] Create `CHANGELOG.md` with template structure

### 1.10 Version Control
- [ ] Commit structure migration: "refactor: migrate to flat package structure"
- [ ] Commit build system: "build: migrate to pyproject.toml"
- [ ] Commit dev infrastructure: "build: add Makefile and update CI"
- [ ] Commit test structure: "test: reorganize test suite"
- [ ] Commit ruff fixes: "style: apply ruff formatting and linting"
- [ ] Commit mypy fixes: "refactor: add type hints for mypy compliance"
- [ ] Commit docstring migration: "docs: migrate to Google-style docstrings"
- [ ] Commit cleanup: "chore: remove deprecated files and update version"

---

## Phase 2: Core Features

### 2.1 Matrix Serialization (YAML)
- [ ] Implement `ScoringMatrix.from_yaml()`
- [ ] Implement `ScoringMatrix.to_yaml()`
- [ ] Implement `ScoringMatrix.from_file()` (auto-detect format)
- [ ] Add YAML dependency to pyproject.toml
- [ ] Update `save()` method to support format parameter
- [ ] Add tests for YAML serialization/deserialization
- [ ] Create example YAML matrices in `tests/data/matrices/`
- [ ] Update documentation with YAML format

### 2.2 Matrix Builder Pattern
- [ ] Refactor existing constructor to `from_scores()`
- [ ] Implement `ScoringMatrix.from_sequences()` (identity matrix)
- [ ] Implement `ScoringMatrix.from_json()` (explicit legacy)
- [ ] Update `__init__()` to be simpler/internal
- [ ] Add `unobserved_strategy` parameter
- [ ] Implement unobserved handling: "min", "impute", "error"
- [ ] Update all tests to use builder methods
- [ ] Update examples to use builder methods

### 2.3 Alignment Validation Metrics
- [ ] Implement `Alignment.sum_of_pairs_score(matrix)`
- [ ] Implement `Alignment.entropy()`
- [ ] Implement `Alignment.consensus(threshold)`
- [ ] Implement `Alignment.gap_statistics()`
- [ ] Implement `malign.sum_of_pairs(alm, matrix)` (standalone)
- [ ] Implement `malign.alignment_entropy(alm)` (standalone)
- [ ] Implement `malign.alignment_consensus(alm, threshold)` (standalone)
- [ ] Implement `malign.alignment_accuracy(alm, gold)` (standalone)
- [ ] Add tests for all metrics
- [ ] Add examples to tutorials

### 2.4 Main API Updates
- [ ] Rename `multi_align()` → `align()`
- [ ] Update `align()` to always return `List[Alignment]`
- [ ] Implement `align_batch()` for multiple alignments
- [ ] Update `tabulate_alms()` for new structure
- [ ] Update `__init__.py` exports
- [ ] Add deprecation warning for old function names (optional)
- [ ] Update all tests for new API
- [ ] Update all examples for new API

### 2.5 Matrix Learning - EM Method
- [ ] Create `malign/matrix_learning.py` module
- [ ] Implement `learn_matrix()` function
- [ ] Implement EM algorithm core loop
- [ ] Implement E-step: align with current matrix
- [ ] Implement M-step: update scores from alignments
- [ ] Implement convergence checking
- [ ] Add support for initial matrix
- [ ] Add support for priors (freqprob integration)
- [ ] Add tests for EM learning
- [ ] Create tutorial example

### 2.6 Matrix Learning - Gradient Descent (if time permits)
- [ ] Implement gradient descent optimizer
- [ ] Define differentiable loss function
- [ ] Implement gradient computation
- [ ] Implement parameter updates with constraints
- [ ] Add learning rate scheduling
- [ ] Add tests for gradient descent
- [ ] Compare with EM in benchmarks

### 2.7 Cross-Validation
- [ ] Implement `cross_validate_matrix()`
- [ ] Implement `split_cognates()` helper
- [ ] Implement `evaluate_matrix()` helper
- [ ] Add tests for CV functions
- [ ] Add tutorial example

### 2.8 External Integrations
- [ ] Implement `ScoringMatrix.from_freqprob()`
- [ ] Implement `ScoringMatrix.from_asymcat()`
- [ ] Add freqprob to optional dependencies
- [ ] Add asymcat to optional dependencies
- [ ] Add integration tests (if libraries available)
- [ ] Add examples to tutorials

---

## Phase 3: Testing & Validation

### 3.1 Test Coverage
- [ ] Run `make test-cov` and assess current coverage
- [ ] Identify uncovered code paths
- [ ] Write tests to reach >90% coverage
- [ ] Add edge case tests (empty, single, all-gap sequences)
- [ ] Add error handling tests (invalid inputs)
- [ ] Verify coverage threshold in CI

### 3.2 Integration Tests
- [ ] Test: sequences → identity matrix → align
- [ ] Test: sequences → learned matrix → align
- [ ] Test: align → validation metrics
- [ ] Test: align_batch with multiple matrices
- [ ] Test: YAML save → load → align (round-trip)
- [ ] Test: freqprob integration (if available)
- [ ] Test: asymcat integration (if available)

### 3.3 Property-Based Tests
- [ ] Test: alignment length consistency
- [ ] Test: symbol preservation (no gaps)
- [ ] Test: score ordering (best first)
- [ ] Test: k-best uniqueness
- [ ] Test: serialization round-trip (matrix)
- [ ] Test: commutativity (where applicable)

### 3.4 Regression Tests (requires gold data)
- [ ] Obtain gold-standard cognate sets from maintainer
- [ ] Format gold data as YAML
- [ ] Implement regression test suite
- [ ] Test against known alignments
- [ ] Set accuracy thresholds (e.g., 85%)
- [ ] Add to CI pipeline

### 3.5 Performance Benchmarks
- [ ] Create `scripts/benchmarks.py`
- [ ] Benchmark ANW vs YenKSP (various k)
- [ ] Benchmark small vs medium workloads
- [ ] Benchmark matrix learning convergence
- [ ] Determine computational limits (max sequences)
- [ ] Document results in USER_GUIDE.md
- [ ] Create performance comparison plots

### 3.6 Comparative Tests (vs lingpy, if possible)
- [ ] Install lingpy (if available)
- [ ] Identify comparable functionality
- [ ] Implement side-by-side comparison
- [ ] Document where MAlign excels
- [ ] Create comparison tutorial/example

---

## Phase 4: Documentation

### 4.1 Nhandu Tutorials
- [ ] Complete `tutorial_1_basics.py`
- [ ] Complete `tutorial_2_matrices.py`
- [ ] Complete `tutorial_3_learning.py`
- [ ] Complete `tutorial_4_advanced.py`
- [ ] Generate HTML: `make docs`
- [ ] Verify all outputs render correctly
- [ ] Add figures to `docs/figures/`

### 4.2 API Reference
- [ ] Auto-generate from docstrings (pydoc-markdown or similar)
- [ ] Review and edit `docs/API_REFERENCE.md`
- [ ] Ensure all public functions documented
- [ ] Add cross-references between functions
- [ ] Add usage examples for each function

### 4.3 User Guide
- [ ] Write comprehensive `docs/USER_GUIDE.md`
- [ ] Section 1: Getting Started
- [ ] Section 2: Core Concepts
- [ ] Section 3: Matrix Management
- [ ] Section 4: Algorithms (ANW vs YenKSP guidance)
- [ ] Section 5: Advanced Topics
- [ ] Section 6: Performance & Limits
- [ ] Section 7: Troubleshooting
- [ ] Add benchmarking results
- [ ] Add decision flowcharts

### 4.4 LLM Documentation
- [ ] Create single-file `docs/LLM_DOCUMENTATION.md`
- [ ] Include all API specifications
- [ ] Include design philosophy
- [ ] Include complete examples
- [ ] Include common patterns
- [ ] Test with LLM (verify comprehension)

### 4.5 README Update
- [ ] Update README.md for v0.4.0
- [ ] Highlight breaking changes
- [ ] Update installation instructions
- [ ] Update quick start examples
- [ ] Add links to tutorials
- [ ] Add performance notes
- [ ] Add badges (CI, PyPI, coverage)

### 4.6 Example Gallery
- [ ] DNA alignment example (update existing)
- [ ] Cross-linguistic example (Italian-Russian)
- [ ] Cognate set learning example
- [ ] Asymmetry showcase example
- [ ] Batch processing example
- [ ] Custom matrix construction example

---

## Phase 5: Release Preparation

### 5.1 Changelog
- [ ] Create `CHANGELOG.md` following Keep a Changelog
- [ ] Document all breaking changes
- [ ] Document all new features
- [ ] Document bug fixes
- [ ] Add migration guide (0.3.x → 0.4.0)

### 5.2 Version Management
- [ ] Update `__version__` to "0.4.0"
- [ ] Update `pyproject.toml` version (if separate)
- [ ] Tag version in git: `git tag v0.4.0`

### 5.3 Quality Assurance
- [ ] Run full test suite: `make test-cov`
- [ ] Verify >90% coverage
- [ ] Run quality checks: `make quality`
- [ ] Verify ruff passes
- [ ] Verify mypy passes
- [ ] Test in clean virtual environment
- [ ] Test installation: `pip install -e .`
- [ ] Test build: `make build`

### 5.4 Documentation Review
- [ ] Proofread all documentation
- [ ] Verify all links work
- [ ] Verify all examples run
- [ ] Verify nhandu tutorials generate correctly
- [ ] Check for typos and clarity

### 5.5 Pre-Release Testing
- [ ] Share with 2+ beta testers
- [ ] Collect feedback on API usability
- [ ] Address critical issues
- [ ] Update documentation based on feedback

---

## Phase 6: Publication

### 6.1 PyPI Release
- [ ] Build distribution: `make build-release`
- [ ] Test upload to TestPyPI: `twine upload --repository testpypi dist/*`
- [ ] Verify TestPyPI installation
- [ ] Upload to PyPI: `twine upload dist/*`
- [ ] Verify PyPI installation: `pip install malign==0.4.0`

### 6.2 GitHub Release
- [ ] Push all changes: `git push`
- [ ] Push tags: `git push --tags`
- [ ] Create GitHub release for v0.4.0
- [ ] Upload built distributions as artifacts
- [ ] Copy CHANGELOG.md entry to release notes

### 6.3 Announcement
- [ ] Announce on relevant mailing lists (linguistics)
- [ ] Tweet/social media (if applicable)
- [ ] Update project homepage (if exists)
- [ ] Notify beta testers

### 6.4 JOSS Preparation (Future)
- [ ] Review JOSS submission guidelines
- [ ] Prepare paper draft
- [ ] Conduct comparison study vs lingpy
- [ ] Collect validation results
- [ ] Submit to JOSS (separate effort, post-release)

---

## Success Metrics

### Code Quality (Phase 1-3)
- ✓ 100% pass on ruff formatting and linting
- ✓ 100% pass on mypy type checking
- ✓ >90% test coverage (enforced)
- ✓ Zero known critical bugs

### Functionality (Phase 2-3)
- ✓ Matrix learning (EM method) functional
- ✓ YAML serialization working
- ✓ Batch API functional
- ✓ Validation metrics implemented
- ✓ Builder pattern for matrices
- ✓ Integration with freqprob/asymcat (basic)

### Documentation (Phase 4)
- ✓ 4+ complete nhandu tutorials
- ✓ Comprehensive USER_GUIDE.md
- ✓ Complete API_REFERENCE.md
- ✓ 6+ working examples

### Performance (Phase 3)
- ✓ Computational limits documented
- ✓ ANW vs YenKSP benchmarked
- ✓ No performance regressions vs v0.3.x

### Community (Phase 6)
- ✓ PyPI release successful
- ✓ 2+ beta testers validated pre-release
- ✓ Positive feedback on API clarity

---

## Timeline Estimates

**Phase 1** (Structure & Quality): 2-3 weeks
- Structure migration: 3-4 days
- Ruff compliance: 3-4 days
- MyPy compliance: 5-7 days
- Docstring migration: 3-4 days
- Buffer: 2-3 days

**Phase 2** (Core Features): 4-6 weeks
- YAML + builders: 1 week
- Validation metrics: 1 week
- Matrix learning (EM): 2-3 weeks
- Integrations: 1 week
- Buffer: 1 week

**Phase 3** (Testing): 2-3 weeks
- Coverage: 1 week
- Property/integration tests: 1 week
- Benchmarks: 3-5 days
- Buffer: 2-3 days

**Phase 4** (Documentation): 3-4 weeks
- Tutorials: 2 weeks
- API reference: 3-4 days
- User guide: 1 week
- Buffer: 3-4 days

**Phase 5** (Release Prep): 1 week
- Changelog: 1 day
- QA: 2-3 days
- Beta testing: 2-3 days

**Phase 6** (Publication): 1 week
- PyPI release: 1 day
- Announcements: 1-2 days
- Buffer: 3-4 days

**Total**: 13-18 weeks focused work (3-4.5 months)
**Calendar**: 6-12 months for side project

**Critical Path**: Matrix learning implementation + validation

---

## Next Steps

1. ✓ **Review and approve this specification** (DONE)
2. **Begin Phase 1**: Code quality & structure migration
   - Start with: Project structure migration (1.1)
   - Then: Build system migration (1.2)
3. **Track progress**: Use TODO list tool for phase tracking
4. **Commit frequently**: Small, focused commits with clear messages

---

## Open Questions (Deferred)

These questions don't block Phase 1, but need answers before later phases:

### For Phase 2:
- **freqprob integration**: Exact API for `from_freqprob()` depends on freqprob internals
- **asymcat integration**: Same for `from_asymcat()`
- **Gradient descent**: Mathematical details of loss function and constraints

### For Phase 3:
- **Computational limits**: Need empirical testing to determine (5? 10? 20 sequences?)
- **ANW vs YenKSP guidance**: Depends on benchmark results

### For Phase 4:
- **Comparison to lingpy**: Depends on comparable functionality and available test data
- **nhandu examples**: May need adjustments based on nhandu updates

### For Phase 5:
- **Beta testers**: Who will test the pre-release?

---

**Document Version**: 2.0
**Last Updated**: 2025-10-14
**Status**: Approved, Ready for Phase 1
