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
- [ ] Add `from_distfeat()` support for directed distances (if distfeat
  provides them), or document the symmetry limitation.

### 2. Different Domains -- DELIVERED

Sequences can have completely different alphabets (Latin, Cyrillic, IPA, etc.).
Each domain is independent. Well-tested and documented.

### 3. True Multi-alignment -- DELIVERED

**What works:**
- For N <= 4 sequences (within grid-size limits), true N-dimensional alignment
  is performed via YenKSP on an N-dimensional graph. All sequences are scored
  jointly at every column, finding the globally optimal alignment.
- For N > 4 or very long sequences, the system silently falls back to the
  progressive pairwise approach for practical feasibility.
- Pairwise alignment (N=2) is dispatched directly without overhead.

**Architecture:**
- `ndim_common.py`: Move generation and feasibility thresholds.
- `ndim_yenksp.py`: N-dimensional graph construction and YenKSP alignment.
- `ndim_nw.py`: N-dimensional Needleman-Wunsch (available but not used in
  dispatch due to exponential co-optimal path explosion with identity matrices).

**Remaining work:**
- [ ] Consider UPGMA-guided progressive alignment for large N (better than
  current all-pairs progressive approach).
- [ ] Investigate capping NW backtrace to enable N-dim NW for cases where
  scoring matrices have few ties.

### 4. Automatic Matrix Inference -- PARTIAL

**What works:**
- `learn_matrix()` learns from pre-clustered cognate sets (EM or gradient
  descent).
- `from_distfeat()` infers scores from phonological feature knowledge.
- `from_substitution_counts()` converts observed counts to log-odds scores.
- Matrix imputation fills sparse matrices from partial data.

**What's missing:**
- No unsupervised mode: users must provide cognate groupings.
- No automatic clustering or cognate discovery.

**Planned:**
- [ ] Consider adding an unsupervised bootstrap that treats each sequence pair
  as its own cognate set, iteratively aligns and re-estimates.
- [ ] Document the supervised requirement clearly.
