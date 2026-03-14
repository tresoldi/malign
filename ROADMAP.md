# MAlign Roadmap

Status of key selling points and planned work.

## Selling Points

### 1. Asymmetric Scoring

- The data structure stores directional scores as ordered tuples
  (`("A","B")` can differ from `("B","A")`). Both ANW and YenKSP use directional
  lookups.
- `ScoringMatrix.from_substitution_counts()` creates asymmetric matrices from
  observed substitution frequencies using log-odds scoring (same approach as
  BLOSUM/PAM matrices).
- N-dimensional alignment natively handles asymmetric matrices.
- `from_distfeat()` produces symmetric feature distances; for asymmetric scoring,
  use `from_substitution_counts()` with directional substitution data, or refine
  a distfeat prior via `bootstrap_matrix()`/`learn_matrix()`.

### 2. Different Domains

Sequences can have completely different alphabets (Latin, Cyrillic, IPA, etc.).
Each domain is independent. Well-tested and documented.

### 3. True Multi-alignment

- For N <= 4 sequences (within grid-size limits), true N-dimensional alignment
  is performed via YenKSP on an N-dimensional graph.
- For N > 4 or very long sequences, UPGMA-guided progressive alignment with
  beam search for k-best results.
- Pairwise alignment (N=2) is dispatched directly without overhead.

### 4. Automatic Matrix Inference

- `learn_matrix()` learns from pre-clustered cognate sets (EM or gradient
  descent). Supports `prior_matrix` for regularization with linear decay.
- `bootstrap_matrix()` learns from arbitrary sequence pairs without clustering
  (unsupervised). Uses iterative log-odds re-estimation with Laplace smoothing.
- Both `learn_matrix()` and `bootstrap_matrix()` accept a `prior_matrix`
  parameter (e.g. from `from_distfeat()`) whose influence decays linearly
  over iterations, enabling the **distfeat → learn/bootstrap → asymmetric**
  pipeline.
- `from_distfeat()` infers scores from phonological feature knowledge
  (symmetric only).
- `from_substitution_counts()` converts observed counts to log-odds scores
  (supports asymmetric).
- Matrix imputation fills sparse matrices from partial data.

### 5. Block Detection and Column Merging

- `detect_blocks()` identifies complementary-gap patterns (diphthongization,
  metathesis) in alignments. Configurable `max_block_size` (default 2).
- `merge_alignment_blocks()` merges block columns into compound symbols
  (tuples). Available as post-processing via `align(merge_blocks=True)`.
- `bootstrap_matrix(block_merge=True)` adjusts pair counts to reduce gap
  inflation caused by block patterns during matrix learning.

## Open Work

- [ ] Investigate capping NW backtrace to enable N-dim NW for cases where
  scoring matrices have few ties.
- [ ] Automatic clustering or cognate discovery (beyond current scope).
