# MAlign Roadmap

Status of key selling points and planned work to close gaps.

## Selling Points Assessment

### 1. Asymmetric Scoring -- DELIVERED

**What works:**
- The data structure stores directional scores as ordered tuples
  (`("A","B")` can differ from `("B","A")`). Both ANW and YenKSP use directional
  lookups.
- `ScoringMatrix.from_substitution_counts()` creates asymmetric matrices from
  observed substitution frequencies using log-odds scoring (same approach as
  BLOSUM/PAM matrices).
- N-dimensional alignment natively handles asymmetric matrices.

**Remaining work:**
- [x] Document `from_distfeat()` symmetry limitation (distfeat provides
  symmetric feature distances; for asymmetric scoring, use
  `from_substitution_counts()` with directional substitution data).

### 2. Different Domains -- DELIVERED

Sequences can have completely different alphabets (Latin, Cyrillic, IPA, etc.).
Each domain is independent. Well-tested and documented.

### 3. True Multi-alignment -- DELIVERED

**What works:**
- For N <= 4 sequences (within grid-size limits), true N-dimensional alignment
  is performed via YenKSP on an N-dimensional graph. All sequences are scored
  jointly at every column, finding the globally optimal alignment.
- For N > 4 or very long sequences, the system uses UPGMA-guided progressive
  alignment with beam search for k-best results.
- Pairwise alignment (N=2) is dispatched directly without overhead.

**Architecture:**
- `ndim_common.py`: Move generation and feasibility thresholds.
- `ndim_yenksp.py`: N-dimensional graph construction and YenKSP alignment.
- `ndim_nw.py`: N-dimensional Needleman-Wunsch (available but not used in
  dispatch due to exponential co-optimal path explosion with identity matrices).
- `progressive.py`: UPGMA-guided progressive alignment with beam search.

**Remaining work:**
- [x] UPGMA-guided progressive alignment for large N (replaces all-pairs
  progressive approach).
- [ ] Investigate capping NW backtrace to enable N-dim NW for cases where
  scoring matrices have few ties.

### 4. Automatic Matrix Inference -- DELIVERED

**What works:**
- `learn_matrix()` learns from pre-clustered cognate sets (EM or gradient
  descent). This is supervised and requires pre-grouped cognates.
- `bootstrap_matrix()` learns from arbitrary sequence pairs without clustering
  (unsupervised). Uses iterative log-odds re-estimation with Laplace smoothing.
- `from_distfeat()` infers scores from phonological feature knowledge
  (symmetric only).
- `from_substitution_counts()` converts observed counts to log-odds scores
  (supports asymmetric).
- Matrix imputation fills sparse matrices from partial data.

- `bootstrap_matrix(prior_matrix=...)` blends a phonological feature prior
  (e.g. from `from_distfeat()`) with data-driven log-odds throughout
  learning. The prior weight decays linearly so phonological structure
  dominates early and fades as the data signal strengthens. This enables
  the key workflow: **distfeat → bootstrap → asymmetric**, starting from
  universal phonological knowledge and refining with language-specific pairs.

**Remaining work:**
- [ ] Automatic clustering or cognate discovery (beyond current scope).

### 5. Block Detection and Column Merging -- DELIVERED

**What works:**
- `detect_blocks()` identifies complementary-gap patterns (diphthongization,
  metathesis) in alignments. Configurable `max_block_size` (default 2).
- `merge_alignment_blocks()` merges block columns into compound symbols
  (tuples). Available as post-processing via `align(merge_blocks=True)`.
- `bootstrap_matrix(block_merge=True)` adjusts pair counts to reduce gap
  inflation caused by block patterns during matrix learning.
- Broadened block definition: single PARTIAL columns extend to adjacent
  FULL columns to capture pairwise diphthongization patterns.
