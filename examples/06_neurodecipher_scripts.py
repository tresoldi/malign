"""Cross-script alignment: Greek -> Linear B, with bootstrap learning.

Demonstrates: align(), ScoringMatrix.from_sequences(), bootstrap_matrix().

Data: data/neurodecipher/greek_linearb.tsv (Greek words paired with
Linear B syllabary transcriptions).
"""

import csv
import random
from pathlib import Path

from malign import align, bootstrap_matrix

DATA_DIR = Path(__file__).parent / "data"
OUTPUT_DIR = Path(__file__).parent / "output"
OUTPUT_DIR.mkdir(exist_ok=True)

# --- load data ---
pairs: list[tuple[list[str], list[str]]] = []
with open(DATA_DIR / "neurodecipher" / "greek_linearb.tsv", encoding="utf-8") as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    for row in reader:
        greek = row["GREEK"].strip()
        linearb = row["LINEAR_B"].strip()
        if greek and linearb:
            # Split words into individual characters
            greek_chars = list(greek)
            linearb_chars = list(linearb)
            pairs.append((greek_chars, linearb_chars))

print(f"Loaded {len(pairs)} Greek-Linear B pairs")

# Sample for speed
random.seed(42)
sample = random.sample(pairs, min(200, len(pairs)))

# --- Step 1: Align with identity matrix ---
lines: list[str] = ["=== Greek -> Linear B cross-script alignment ===", ""]
lines.append("--- Identity matrix (first 10 pairs) ---")
for greek, linearb in sample[:10]:
    alms = align([greek, linearb], k=1)
    if alms:
        a = alms[0]
        lines.append(f"  Greek:    {''.join(str(s) for s in a.seqs[0])}")
        lines.append(f"  Linear B: {''.join(str(s) for s in a.seqs[1])}")
        lines.append(f"  score: {a.score:.2f}")
        lines.append("")

# --- Step 2: Bootstrap a scoring matrix ---
print("Bootstrapping scoring matrix from 200 pairs...")
learned = bootstrap_matrix(
    sample,
    max_iter=5,
    verbose=True,
)

# --- Step 3: Re-align with learned matrix ---
lines.append("--- Learned matrix (first 10 pairs) ---")
for greek, linearb in sample[:10]:
    alms = align([greek, linearb], matrix=learned, k=1)
    if alms:
        a = alms[0]
        lines.append(f"  Greek:    {''.join(str(s) for s in a.seqs[0])}")
        lines.append(f"  Linear B: {''.join(str(s) for s in a.seqs[1])}")
        lines.append(f"  score: {a.score:.2f}")
        lines.append("")

alm_file = OUTPUT_DIR / "neurodecipher_alignments.txt"
alm_file.write_text("\n".join(lines), encoding="utf-8")
print(f"Alignments written to {alm_file}")

# --- Write top cross-script scores ---
score_lines: list[str] = ["=== Top 30 Greek-Linear B character scores ===", ""]
# Filter to cross-script pairs (not gap-gap)
gap = learned.gap
cross_scores = [
    (pair, score)
    for pair, score in learned.scores.items()
    if gap not in pair
]
cross_scores.sort(key=lambda x: x[1], reverse=True)
for pair, score in cross_scores[:30]:
    score_lines.append(f"  {pair[0]} -> {pair[1]}  {score:+.3f}")

top_file = OUTPUT_DIR / "neurodecipher_top_scores.txt"
top_file.write_text("\n".join(score_lines), encoding="utf-8")
print(f"Top scores written to {top_file}")
