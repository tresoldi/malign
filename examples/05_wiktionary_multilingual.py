"""Multilingual alignment of Romance cognates from Wiktionary.

Demonstrates: align(), bootstrap_matrix(), comparison of identity vs
learned matrix, and 3-way alignment.

Italian and Portuguese are closely related Romance languages, so most
translation equivalents share common Latin etyma and produce meaningful
alignments. English is added as a third language for 3-way demos
(shared cognates via Latin/French borrowings).

Data: data/wiktionary/wiktionary.tsv.
"""

import csv
import random
from pathlib import Path

from malign import align, bootstrap_matrix

DATA_DIR = Path(__file__).parent / "data"
OUTPUT_DIR = Path(__file__).parent / "output"
OUTPUT_DIR.mkdir(exist_ok=True)

# --- load data ---
pairs_it_pt: list[tuple[list[str], list[str]]] = []
triples: list[tuple[list[str], list[str], list[str]]] = []

with open(DATA_DIR / "wiktionary" / "wiktionary.tsv", encoding="utf-8") as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    for row in reader:
        it = row.get("Italian", "").strip()
        pt = row.get("Portuguese", "").strip()
        if not it or not pt:
            continue
        it_chars = it.split()
        pt_chars = pt.split()
        # Keep short single-word entries (likely cognates, not compounds)
        if len(it_chars) > 12 or len(pt_chars) > 12:
            continue
        pairs_it_pt.append((it_chars, pt_chars))

        # Also collect triples with English
        en = row.get("English", "").strip()
        if en:
            en_chars = en.split()
            if len(en_chars) <= 12:
                triples.append((it_chars, pt_chars, en_chars))

print(f"Loaded {len(pairs_it_pt)} Italian-Portuguese pairs, {len(triples)} triples with English")

# Reproducible sample
random.seed(42)
pair_sample = random.sample(pairs_it_pt, min(200, len(pairs_it_pt)))
triple_sample = random.sample(triples, min(20, len(triples)))

# --- Step 1: Bootstrap matrix from Italian-Portuguese pairs ---
boot_subset = pair_sample[:150]
print(f"Bootstrapping matrix from {len(boot_subset)} Italian-Portuguese pairs...")
learned = bootstrap_matrix(
    boot_subset,
    max_iter=10,
    verbose=True,
)

# --- Step 2: Compare identity vs learned on pairwise alignments ---
lines: list[str] = [
    "=== Italian-Portuguese alignment (identity vs bootstrapped) ===",
    "",
]

for it_chars, pt_chars in pair_sample[:20]:
    it_word = "".join(it_chars)
    pt_word = "".join(pt_chars)
    lines.append(f"--- {it_word} / {pt_word} ---")

    for label, matrix in [("identity", None), ("learned", learned)]:
        alms = align([it_chars, pt_chars], k=1, matrix=matrix)
        if alms:
            a = alms[0]
            lines.append(f"  {label:>10}: {' '.join(str(s) for s in a.seqs[0])}")
            lines.append(f"  {'':>10}  {' '.join(str(s) for s in a.seqs[1])}")
            lines.append(f"  {'':>10}  score: {a.score:.2f}")
    lines.append("")

# --- Step 3: 3-way alignment (Italian / Portuguese / English) ---
lines.append("=== 3-way alignment (Italian / Portuguese / English) ===")
lines.append("")

for it_chars, pt_chars, en_chars in triple_sample:
    alms = align([it_chars, pt_chars, en_chars], k=1, matrix=learned)
    if not alms:
        continue
    a = alms[0]
    lines.append(f"  {'Italian':>12}: {' '.join(str(s) for s in a.seqs[0])}")
    lines.append(f"  {'Portuguese':>12}: {' '.join(str(s) for s in a.seqs[1])}")
    lines.append(f"  {'English':>12}: {' '.join(str(s) for s in a.seqs[2])}")
    lines.append(f"  {'score':>12}: {a.score:.2f}")
    lines.append("")

outfile = OUTPUT_DIR / "wiktionary_multilingual.txt"
outfile.write_text("\n".join(lines), encoding="utf-8")
print(f"\nAlignments written to {outfile}")
