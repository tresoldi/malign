"""Large-scale matrix learning from cross-linguistic data with gold evaluation.

Demonstrates: learn_matrix(), align(), ScoringMatrix.from_sequences(),
alignment_accuracy().

Data: data/lexibank/forms.tsv (filtered to savelyevturkic dataset,
Turkish + Azeri, grouped by cognate class).
"""

import csv
from collections import defaultdict
from pathlib import Path

from malign import Alignment, align, alignment_accuracy, learn_matrix

DATA_DIR = Path(__file__).parent / "data"
OUTPUT_DIR = Path(__file__).parent / "output"
OUTPUT_DIR.mkdir(exist_ok=True)

DATASET = "savelyevturkic"
LANG_A = "savelyevturkic_Turkish"
LANG_B = "savelyevturkic_Azeri"

# --- load data ---
# Group by cognacy: {cognacy_id: {lang_id: (segments, alignment)}}
cognates: dict[str, dict[str, tuple[list[str], list[str]]]] = defaultdict(dict)

with open(DATA_DIR / "lexibank" / "forms.tsv", encoding="utf-8") as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    for row in reader:
        if row["Dataset"] != DATASET:
            continue
        lang = row["Language_ID"]
        if lang not in (LANG_A, LANG_B):
            continue
        cog = row["Cognacy"].strip()
        seg = row["Segments"].strip()
        alm = row["Alignment"].strip()
        if cog and seg:
            cognates[cog][lang] = (seg.split(), alm.split() if alm else [])

# Filter to cognate sets with both languages present
paired: list[tuple[list[str], list[str], list[str], list[str]]] = []
for cog_id, langs in cognates.items():
    if LANG_A in langs and LANG_B in langs:
        seg_a, alm_a = langs[LANG_A]
        seg_b, alm_b = langs[LANG_B]
        paired.append((seg_a, seg_b, alm_a, alm_b))

# Limit to 150 pairs for speed
paired = paired[:150]
print(f"Found {len(paired)} {LANG_A.split('_')[1]}-{LANG_B.split('_')[1]} cognate pairs")

if not paired:
    print("No paired cognates found. Try different languages.")
    raise SystemExit(1)

# --- Step 1: Learn matrix from cognate sets ---
cognate_sets = [[sa, sb] for sa, sb, _, _ in paired]

print("Learning scoring matrix...")
learned = learn_matrix(
    cognate_sets,
    max_iter=5,
    verbose=True,
)

mat_file = OUTPUT_DIR / "lexibank_learned_matrix.txt"
mat_file.write_text(learned.tabulate(), encoding="utf-8")
print(f"Learned matrix written to {mat_file}")

# --- Step 2: Align with learned matrix and evaluate ---
lines: list[str] = [f"=== {LANG_A.split('_')[1]}-{LANG_B.split('_')[1]} alignments ===", ""]
evaluated = 0
total_acc = 0.0

for sa, sb, gold_a, gold_b in paired[:50]:
    alms = align([sa, sb], matrix=learned, k=1)
    if not alms:
        continue

    pred = alms[0]
    lines.append(f"  {' '.join(str(s) for s in pred.seqs[0])}")
    lines.append(f"  {' '.join(str(s) for s in pred.seqs[1])}")
    lines.append(f"  score: {pred.score:.2f}")

    # Evaluate against gold if available and lengths match
    if gold_a and gold_b and len(gold_a) == len(gold_b):
        gold = Alignment(seqs=[gold_a, gold_b], score=0.0)
        if len(pred.seqs[0]) == len(gold.seqs[0]):
            acc = alignment_accuracy(pred, gold)
            total_acc += acc
            evaluated += 1
            lines.append(f"  accuracy vs gold: {acc:.4f}")

    lines.append("")

alm_file = OUTPUT_DIR / "lexibank_alignments.txt"
alm_file.write_text("\n".join(lines), encoding="utf-8")
print(f"Alignments written to {alm_file}")

if evaluated > 0:
    print(f"Evaluated {evaluated} alignments against gold, avg accuracy: {total_acc / evaluated:.4f}")
else:
    print("No gold alignments matched for evaluation.")
