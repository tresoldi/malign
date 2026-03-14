"""Distfeat -> bootstrap_matrix(prior) -> learn_matrix(prior) pipeline.

Demonstrates: ScoringMatrix.from_distfeat(), bootstrap_matrix(),
learn_matrix(), align().

Data: data/northeuralex/italic_by_cogid.tsv (Romance cognates).
Requires: pip install malign[features]
"""

import csv
from pathlib import Path

from malign import align, bootstrap_matrix, learn_matrix
from malign.scoring_matrix import ScoringMatrix

DATA_DIR = Path(__file__).parent / "data"
OUTPUT_DIR = Path(__file__).parent / "output"
OUTPUT_DIR.mkdir(exist_ok=True)

# --- check distfeat availability ---
try:
    ScoringMatrix.from_distfeat([["a", "b"], ["a", "b"]])
    HAS_DISTFEAT = True
except ImportError:
    HAS_DISTFEAT = False
    print("distfeat not installed. Install with: pip install malign[features]")
    print("Falling back to identity matrix as prior.\n")

# --- load data ---
LANG_A, LANG_B = "ITALIAN", "SPANISH"

cognate_pairs: list[tuple[list[str], list[str]]] = []
cognate_sets: list[list[list[str]]] = []

with open(DATA_DIR / "northeuralex" / "italic_by_cogid.tsv", encoding="utf-8") as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    for row in reader:
        seq_a = row.get(f"{LANG_A}_SEQ", "").strip()
        seq_b = row.get(f"{LANG_B}_SEQ", "").strip()
        if seq_a and seq_b:
            # Filter out morpheme boundary markers
            sa = [s for s in seq_a.split() if s != "+"]
            sb = [s for s in seq_b.split() if s != "+"]
            if sa and sb:
                cognate_pairs.append((sa, sb))
                cognate_sets.append([sa, sb])

# Limit to first 100 pairs for speed
cognate_pairs = cognate_pairs[:100]
cognate_sets = cognate_sets[:100]

print(f"Loaded {len(cognate_pairs)} {LANG_A}-{LANG_B} cognate pairs")

# --- Step 1: Build phonological prior (or identity fallback) ---
all_symbols_a: set[str] = set()
all_symbols_b: set[str] = set()
for sa, sb in cognate_pairs:
    all_symbols_a.update(sa)
    all_symbols_b.update(sb)

if HAS_DISTFEAT:
    prior = ScoringMatrix.from_distfeat(
        [sorted(all_symbols_a), sorted(all_symbols_b)],
        impute_method="mean",
    )
    print(f"Built distfeat prior: {len(prior.scores)} score entries")
else:
    prior = ScoringMatrix.from_sequences(
        [sorted(all_symbols_a), sorted(all_symbols_b)],
        match=1.0,
        mismatch=-0.5,
        gap_score=-1.0,
    )
    print(f"Built identity prior: {len(prior.scores)} score entries")

# Save prior matrix
prior_file = OUTPUT_DIR / "northeuralex_distfeat_matrix.txt"
prior_file.write_text(prior.tabulate(), encoding="utf-8")
print(f"Prior matrix written to {prior_file}")

# --- Step 2: Bootstrap matrix from pairs ---
print("\nBootstrapping matrix from cognate pairs...")
bootstrapped = bootstrap_matrix(
    cognate_pairs,
    max_iter=10,
    prior_matrix=prior,
    prior_weight=0.3,
    verbose=True,
)

bootstrapped_file = OUTPUT_DIR / "northeuralex_learned_matrix.txt"
bootstrapped_file.write_text(bootstrapped.tabulate(), encoding="utf-8")
print(f"Bootstrapped matrix written to {bootstrapped_file}")

# --- Step 3: Learn matrix from cognate sets (using bootstrapped as prior) ---
print("\nLearning matrix from cognate sets...")
learned = learn_matrix(
    cognate_sets,
    max_iter=5,
    prior_matrix=bootstrapped,
    prior_weight=0.5,
    verbose=True,
)

# --- Step 4: Align with each matrix stage and compare ---
lines: list[str] = [
    f"=== {LANG_A}-{LANG_B} alignment comparison ===",
    "",
    "Each cognate pair is aligned with:",
    "  (a) identity/distfeat prior",
    "  (b) bootstrapped matrix",
    "  (c) EM-learned matrix",
    "",
]

for sa, sb in cognate_pairs[:20]:
    word_a = "".join(sa)
    word_b = "".join(sb)
    lines.append(f"--- {word_a} / {word_b} ---")

    for label, matrix in [("prior", prior), ("bootstrap", bootstrapped), ("learned", learned)]:
        alms = align([sa, sb], matrix=matrix, k=1)
        if alms:
            a = alms[0]
            lines.append(f"  {label:>10}: {' '.join(str(s) for s in a.seqs[0])}")
            lines.append(f"  {'':>10}  {' '.join(str(s) for s in a.seqs[1])}")
            lines.append(f"  {'':>10}  score: {a.score:.2f}")
    lines.append("")

alm_file = OUTPUT_DIR / "northeuralex_alignments.txt"
alm_file.write_text("\n".join(lines), encoding="utf-8")
print(f"\nSample alignments written to {alm_file}")
