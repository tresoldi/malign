# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.5.0] - 2025

### Added

#### True N-Dimensional Multi-Alignment
- **N-dimensional YenKSP**: For up to 4 sequences (within grid-size limits),
  performs true N-dimensional alignment on an N-dimensional graph, finding the
  globally optimal alignment across all sequences simultaneously.
- **Automatic dispatch**: `align()` transparently selects the best strategy
  (direct pairwise, N-dim YenKSP, or progressive) based on sequence count
  and grid size.

#### UPGMA Progressive Alignment
- **`upgma_progressive_align()`**: Guide-tree-based progressive alignment for
  5+ sequences, using UPGMA clustering on pairwise alignment distances.
- **Beam search**: Maintains top-k candidates during profile merging for
  k-best progressive alignments.

#### Unsupervised Bootstrap Matrix Learning
- **`bootstrap_matrix()`**: Learn scoring matrices from arbitrary sequence
  pairs without pre-clustered cognate sets. Uses iterative log-odds
  re-estimation with Laplace smoothing.
- **Prior-guided learning**: `prior_matrix` parameter blends phonological
  feature priors (e.g. from `from_distfeat()`) with data-driven scores.
  Prior weight decays linearly to zero over iterations.
- **Block-aware learning**: `block_merge=True` detects complementary-gap
  patterns and reduces gap inflation in pair counts during the M-step.

#### Block Detection and Column Merging
- **`detect_blocks()`**: Identifies complementary-gap patterns
  (diphthongization, metathesis) in alignments. Configurable `max_block_size`.
- **`merge_alignment_blocks()`**: Merges block columns into compound symbols
  (tuples). Exported in public API.
- **`align(merge_blocks=True)`**: Post-processing integration for automatic
  block merging after alignment.

#### Feature-Based Scoring
- **`ScoringMatrix.from_distfeat()`**: Builds scoring matrices from
  phonological feature distances using the `distfeat` library.
  Requires `pip install malign[features]`.
- **`ScoringMatrix.from_substitution_counts()`**: Creates asymmetric scoring
  matrices from observed substitution frequencies using log-odds scoring
  (same approach as BLOSUM/PAM matrices).

### Changed
- **Breaking**: Python >= 3.12 required (was 3.10+)
- **Breaking**: `ScoringMatrix` is now a frozen dataclass (immutable)
  - Removed `__setitem__`, `load()`, `copy()` methods
  - Use factory classmethods: `from_yaml()`, `from_sequences()`, `from_distfeat()`
  - `domains` is now `tuple[tuple[...], ...]` instead of `list[list[...]]`
- **Breaking**: Dropped `networkx` dependency -- replaced with pure Python
  Dijkstra + Yen's KSP implementation
- EM learning creates new `ScoringMatrix` per iteration (functional style)
- Replaced `pairwise_iter` with `itertools.pairwise` (Python 3.12+)
- Google-style docstrings throughout
- `Alignment` added to public API exports
- Comprehensive ruff linting/formatting compliance

### Removed
- `networkx` dependency
- `ScoringMatrix.__setitem__()`, `.load()`, `.copy()`
- `pairwise_iter` utility function

---

## [0.4.0-beta.1] - 2025-10-16

### Added

#### Matrix Learning (Phase 2)
- **EM Algorithm**: Expectation-Maximization for learning scoring matrices from cognate sets
- **Gradient Descent**: L-BFGS-B optimization via `scipy.optimize.minimize`
- **Convergence Detection**:
  - Score-based convergence (relative change < 0.001)
  - Matrix-based convergence (Frobenius norm < 0.01)
  - OR logic: stops when either criterion is met
- **Early Stopping**: Patience-based stopping (default: 5 iterations without improvement)
- **Parameter Bounds**: Hard constraints [-10, 10] for gradient descent
- **Verbose Logging**: Optional convergence tracking with `verbose=True`

#### Serialization (Phase 2)
- **YAML Format**: Primary format for matrix serialization via `save()` and `from_yaml()`
  - Human-readable format with Unicode support
  - Handles asymmetric matrices and multiple domains
  - Preserves gap symbols and domain information
- **Builder Methods**: `ScoringMatrix.from_sequences()` for quick matrix construction

#### Metrics (Phase 2)
- **Alignment Accuracy**: Column-wise matching percentage
- **Precision/Recall**: Position-based evaluation metrics
- **F1 Score**: Harmonic mean of precision and recall

#### Testing Infrastructure (Phase 3)
- **Coverage Tests**: Expanded from 74% to 77% on core malign.py
- **Property-Based Tests**: 6 core properties using Hypothesis
- **Integration Tests**: 5 end-to-end pipeline scenarios
- **Regression Tests**: Validation against Arca Verborum gold standard (451,935 forms)
- **Benchmark Suite**: Performance analysis across sequence count, length, and k values
- **Test Data**: 420 curated cognate sets from Arca Verborum

#### Documentation
- **Algorithm Selection Guide**: Decision flowchart and performance comparison (ANW vs YenKSP)
- **Test Data Documentation**: Comprehensive `tests/data/README.md`
- **Pytest Markers**: `slow`, `integration`, `regression`, `property`, `benchmark`

### Changed

#### Breaking Changes
- **Function Renamed**: `multi_align()` -> `align()`
- **Python Version**: Now requires Python 3.10+ (was 3.7+)
- **Docstring Style**: Migrated from Epytext to Google-style docstrings

#### Improvements
- **Code Quality**: Comprehensive ruff linting and formatting applied
- **Type Annotations**: Added throughout codebase (mypy compatible)
- **Error Messages**: More descriptive validation errors
- **Performance**: Optimized matrix operations in learning algorithms

### Fixed
- **Sorting Issue**: Fixed None value handling in matrix score key sorting
- **Nested If**: Simplified conditional logic in gradient descent
- **Test Stability**: Fixed flaky integration tests with realistic timeouts

### Deprecated
- JSON matrix format (still supported for backward compatibility, but YAML is preferred)

### Migration Guide from 0.3.x

```python
# Old (0.3.x)
from malign import multi_align
alms = multi_align(sequences, k=3, method="anw")

# New (0.4.0+)
from malign import align
alms = align(sequences, k=3, method="anw")
```

---

## [0.3.0] - 2021

### Added
- Type annotations throughout codebase
- Matrix imputation methods
- Support for any hashable Python objects (not just strings)
- Documentation updates

### Changed
- Code improvements and refactoring
- Preparations for public announcement

---

## [0.2.0] - 2020

### Added
- Asymmetric Needleman-Wunsch (ANW) implementation
- Yen's k-shortest paths (YenKSP) implementation
- ScoringMatrix object for asymmetric scoring
- Reproducible alignment sorting

### Changed
- Major revision of alignment algorithms

---

## [0.1.0] - 2019

### Added
- Initial release for internal testing
- Basic multiple sequence alignment
- Community outreach version
