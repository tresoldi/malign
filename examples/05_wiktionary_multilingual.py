"""Multilingual alignment with mixed scripts (3+ sequences).

Demonstrates: align(), ScoringMatrix.from_sequences(), tabulate_alms().
Shows alignment of words across Latin, Bengali, Georgian, Hebrew, and
other scripts.

Data: data/wiktionary/wiktionary.tsv (filtered to rows with >= 3 languages).
"""

import csv
import random
from pathlib import Path

from malign import align

DATA_DIR = Path(__file__).parent / "data"
OUTPUT_DIR = Path(__file__).parent / "output"
OUTPUT_DIR.mkdir(exist_ok=True)

# --- load data ---
LANGUAGES = [
    "English",
    "Bengali",
    "Georgian",
    "Hebrew",
    "Italian",
    "Portuguese",
    "Russian",
    "Finnish",
    "Czech",
    "Swahili",
    "Afrikaans",
    "Danish",
    "Esperanto",
]

rows_with_multi: list[tuple[list[str], list[list[str]]]] = []

with open(DATA_DIR / "wiktionary" / "wiktionary.tsv", encoding="utf-8") as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    for row in reader:
        present = []
        for lang in LANGUAGES:
            val = row.get(lang, "").strip()
            if val:
                present.append((lang, val.split()))
        # Keep rows with >=3 languages and short sequences (max 12 symbols each)
        if len(present) >= 3 and all(len(p[1]) <= 12 for p in present):
            langs_here = [p[0] for p in present]
            seqs_here = [p[1] for p in present]
            rows_with_multi.append((langs_here, seqs_here))

print(f"Found {len(rows_with_multi)} rows with >= 3 languages")

# Sample 50 rows, preferring mixed-script combinations
random.seed(42)

# Prefer rows that include non-Latin scripts
non_latin_langs = {"Bengali", "Georgian", "Hebrew", "Russian"}
mixed_script = [
    (langs, seqs)
    for langs, seqs in rows_with_multi
    if any(lang in non_latin_langs for lang in langs)
]

if len(mixed_script) >= 20:
    sample = random.sample(mixed_script, 20)
else:
    sample = mixed_script + random.sample(
        [r for r in rows_with_multi if r not in mixed_script],
        min(20 - len(mixed_script), len(rows_with_multi) - len(mixed_script)),
    )

print(f"Selected {len(sample)} entries (preferring mixed-script rows)")

# --- align ---
lines: list[str] = ["=== Multilingual alignment (3+ sequences, mixed scripts) ===", ""]

for i, (langs, seqs) in enumerate(sample):
    # Limit to 3 sequences max for tractability
    if len(seqs) > 3:
        seqs = seqs[:3]
        langs = langs[:3]

    alms = align(seqs, k=1)
    if not alms:
        continue

    a = alms[0]
    lines.append(f"--- Entry {i + 1} ({', '.join(langs)}) ---")
    for lang, seq in zip(langs, a.seqs, strict=False):
        lines.append(f"  {lang:>12}: {' '.join(str(s) for s in seq)}")
    lines.append(f"  {'score':>12}: {a.score:.2f}")
    lines.append("")

outfile = OUTPUT_DIR / "wiktionary_multilingual.txt"
outfile.write_text("\n".join(lines), encoding="utf-8")
print(f"\nAlignments written to {outfile}")
