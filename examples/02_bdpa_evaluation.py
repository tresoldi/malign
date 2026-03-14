"""Gold-standard evaluation of predicted alignments against BDPA data.

Demonstrates: align(), Alignment(), alignment_accuracy(),
alignment_precision_recall(), alignment_f1().

Data: data/bdpa/romance.tsv (31 cognate sets, 8 Romance varieties).
"""

import csv
import itertools
from pathlib import Path

from malign import Alignment, align, alignment_accuracy, alignment_f1, alignment_precision_recall

DATA_DIR = Path(__file__).parent / "data"
OUTPUT_DIR = Path(__file__).parent / "output"
OUTPUT_DIR.mkdir(exist_ok=True)

# --- load BDPA romance data ---
rows: list[dict[str, str]] = []
with open(DATA_DIR / "bdpa" / "romance.tsv", encoding="utf-8") as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    for row in reader:
        rows.append(row)

# Extract language names from columns (every lang has _REF and _SEQ)
header_langs = []
for col in rows[0]:
    if col.endswith("_SEQ"):
        header_langs.append(col.removesuffix("_SEQ"))

results: list[dict[str, str]] = []

for row in rows:
    gloss = row["ID"]

    # Collect languages that have data in this row
    available = []
    for lang in header_langs:
        ref = row.get(f"{lang}_REF", "").strip()
        seq = row.get(f"{lang}_SEQ", "").strip()
        if ref and seq:
            available.append(lang)

    # Pairwise evaluation
    for lang_a, lang_b in itertools.combinations(available, 2):
        seq_a = row[f"{lang_a}_SEQ"].split()
        seq_b = row[f"{lang_b}_SEQ"].split()
        ref_a = row[f"{lang_a}_REF"].split()
        ref_b = row[f"{lang_b}_REF"].split()

        # Predict alignment
        alms = align([seq_a, seq_b], k=1)
        if not alms:
            continue
        predicted = alms[0]

        # Build gold alignment from _REF columns
        gold = Alignment(seqs=[ref_a, ref_b], score=0.0)

        # Skip if lengths don't match (different alignment structure)
        if len(predicted.seqs[0]) != len(gold.seqs[0]):
            continue

        acc = alignment_accuracy(predicted, gold)
        prec, rec = alignment_precision_recall(predicted, gold)
        f1 = alignment_f1(predicted, gold)

        results.append(
            {
                "Gloss": gloss,
                "Lang_A": lang_a,
                "Lang_B": lang_b,
                "Accuracy": f"{acc:.4f}",
                "Precision": f"{prec:.4f}",
                "Recall": f"{rec:.4f}",
                "F1": f"{f1:.4f}",
            }
        )

# Write results
outfile = OUTPUT_DIR / "bdpa_evaluation.tsv"
fieldnames = ["Gloss", "Lang_A", "Lang_B", "Accuracy", "Precision", "Recall", "F1"]
with open(outfile, "w", encoding="utf-8", newline="") as fh:
    writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()
    writer.writerows(results)

# Summary
if results:
    avg_acc = sum(float(r["Accuracy"]) for r in results) / len(results)
    avg_f1 = sum(float(r["F1"]) for r in results) / len(results)
    print(f"Evaluated {len(results)} pairwise alignments across {len(rows)} cognate sets")
    print(f"Average accuracy: {avg_acc:.4f}")
    print(f"Average F1:       {avg_f1:.4f}")
else:
    print("No matching-length alignments found for evaluation.")

print(f"\nResults written to {outfile}")
