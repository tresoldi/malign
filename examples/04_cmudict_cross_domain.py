"""Cross-domain alignment: English letters <-> phonemes.

Demonstrates: align(), ScoringMatrix.from_sequences(),
ScoringMatrix.from_substitution_counts().

Data: data/cmudict/cmudict.tsv (sample of 500 entries).
"""

import csv
import random
from collections import Counter
from pathlib import Path

from malign import align
from malign.scoring_matrix import ScoringMatrix

DATA_DIR = Path(__file__).parent / "data"
OUTPUT_DIR = Path(__file__).parent / "output"
OUTPUT_DIR.mkdir(exist_ok=True)

# --- load data ---
entries: list[tuple[list[str], list[str]]] = []
with open(DATA_DIR / "cmudict" / "cmudict.tsv", encoding="utf-8") as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    for row in reader:
        letters = row["Orthography"].split()
        phones = row["Segments"].split()
        if letters and phones:
            entries.append((letters, phones))

# Reproducible sample of 500 entries
random.seed(42)
sample = random.sample(entries, min(500, len(entries)))
print(f"Sampled {len(sample)} entries from {len(entries)} total")

# --- Step 1: Align with simple identity matrix ---
lines: list[str] = ["=== Cross-domain alignment (letters <-> phonemes) ===", ""]

lines.append("--- Identity matrix (first 10 words) ---")
for letters, phones in sample[:10]:
    alms = align([letters, phones], k=1)
    if alms:
        a = alms[0]
        lines.append(f"  {' '.join(str(s) for s in a.seqs[0])}")
        lines.append(f"  {' '.join(str(s) for s in a.seqs[1])}")
        lines.append(f"  score: {a.score:.2f}")
        lines.append("")

# --- Step 2: Build substitution counts from initial alignments ---
print("Aligning sample to collect substitution counts...")
sub_counts: Counter[tuple[str, str]] = Counter()
for letters, phones in sample:
    alms = align([letters, phones], k=1)
    if alms:
        for l_sym, p_sym in zip(alms[0].seqs[0], alms[0].seqs[1], strict=False):
            sub_counts[(str(l_sym), str(p_sym))] += 1

# Build log-odds matrix from counts
matrix = ScoringMatrix.from_substitution_counts(sub_counts)
print(f"Built substitution-count matrix: {len(matrix.scores)} entries")

# --- Step 3: Re-align with learned matrix ---
lines.append("--- Substitution-count matrix (first 10 words) ---")
for letters, phones in sample[:10]:
    alms = align([letters, phones], matrix=matrix, k=1)
    if alms:
        a = alms[0]
        lines.append(f"  {' '.join(str(s) for s in a.seqs[0])}")
        lines.append(f"  {' '.join(str(s) for s in a.seqs[1])}")
        lines.append(f"  score: {a.score:.2f}")
        lines.append("")

# Write alignments
alm_file = OUTPUT_DIR / "cmudict_alignments.txt"
alm_file.write_text("\n".join(lines), encoding="utf-8")
print(f"Alignments written to {alm_file}")

# Write top substitution scores
mat_lines = ["=== Top 30 letter-phoneme substitution scores ===", ""]
top_pairs = sorted(matrix.scores.items(), key=lambda x: x[1], reverse=True)[:30]
for pair, score in top_pairs:
    mat_lines.append(f"  {pair[0]:>3} -> {pair[1]:<5}  {score:+.3f}")

mat_file = OUTPUT_DIR / "cmudict_matrix_sample.txt"
mat_file.write_text("\n".join(mat_lines), encoding="utf-8")
print(f"Matrix sample written to {mat_file}")
