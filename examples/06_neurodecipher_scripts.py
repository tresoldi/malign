"""Cross-script alignment: Greek alphabet vs Linear B syllabary.

Demonstrates bootstrap learning when the two domains share zero symbols.
Greek is alphabetic (one character per sound) while Linear B is a
syllabary (each sign = CV syllable), so e.g. Greek "agora" (a-g-o-r-a)
maps to Linear B a-go-ra (3 signs). This length mismatch means
alignments will have structural gaps, but bootstrap_matrix() can still
discover which Greek characters tend to co-occur with which Linear B
signs.

API: align(), bootstrap_matrix().
Data: data/neurodecipher/greek_linearb.tsv.
"""

import csv
from collections import Counter
from pathlib import Path

from malign import align, bootstrap_matrix

DATA_DIR = Path(__file__).parent / "data"
OUTPUT_DIR = Path(__file__).parent / "output"
OUTPUT_DIR.mkdir(exist_ok=True)

# --- load data ---
all_pairs: list[tuple[list[str], list[str]]] = []
with open(DATA_DIR / "neurodecipher" / "greek_linearb.tsv", encoding="utf-8") as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    for row in reader:
        greek = row["GREEK"].strip()
        linearb = row["LINEAR_B"].strip()
        if greek and linearb:
            all_pairs.append((list(greek), list(linearb)))

print(f"Loaded {len(all_pairs)} Greek-Linear B pairs")

# --- Step 1: Bootstrap with all data ---
print(f"Bootstrapping scoring matrix from {len(all_pairs)} pairs...")
learned = bootstrap_matrix(
    all_pairs,
    max_iter=15,
    verbose=True,
)

# --- Step 2: Show alignments on short pairs (less length mismatch) ---
short_pairs = [(g, lb) for g, lb in all_pairs if len(g) <= 5 and len(lb) >= 2]

lines: list[str] = [
    "=== Greek (alphabet) -> Linear B (syllabary) alignment ===",
    "",
    "Linear B is a syllabary: each sign encodes a CV syllable.",
    "Greek words are longer (more characters) than their Linear B",
    "equivalents, so gaps are structurally expected.",
    "",
    "--- Short-word alignments (learned matrix) ---",
    "",
]

for greek, linearb in short_pairs[:25]:
    alms = align([greek, linearb], matrix=learned, k=1)
    if not alms:
        continue
    a = alms[0]
    lines.append(f"  Greek:    {' '.join(str(s) for s in a.seqs[0])}")
    lines.append(f"  Linear B: {' '.join(str(s) for s in a.seqs[1])}")
    lines.append(f"  score: {a.score:.2f}")
    lines.append("")

alm_file = OUTPUT_DIR / "neurodecipher_alignments.txt"
alm_file.write_text("\n".join(lines), encoding="utf-8")
print(f"Alignments written to {alm_file}")

# --- Step 3: Learned correspondences ---
# Count symbol frequencies to filter noise from rare symbols
greek_freq: Counter[str] = Counter()
linearb_freq: Counter[str] = Counter()
for greek, linearb in all_pairs:
    greek_freq.update(greek)
    linearb_freq.update(linearb)

# Keep symbols that appear at least 10 times
common_greek = {s for s, c in greek_freq.items() if c >= 10}
common_linearb = {s for s, c in linearb_freq.items() if c >= 10}

gap = learned.gap
cross_scores = [
    (pair, score)
    for pair, score in learned.scores.items()
    if gap not in pair and pair[0] in common_greek and pair[1] in common_linearb
]
cross_scores.sort(key=lambda x: x[1], reverse=True)

# For each Greek character, show its best Linear B match
best_per_greek: dict[str, tuple[str, float]] = {}
for pair, score in cross_scores:
    g_char = pair[0]
    if g_char not in best_per_greek:
        best_per_greek[g_char] = (pair[1], score)

# Known Linear B phonetic values for reference
LB_VALUES = {
    "𐀀": "a", "𐀁": "e", "𐀂": "i", "𐀃": "o", "𐀄": "u",
    "𐀅": "da", "𐀆": "de", "𐀇": "di", "𐀈": "do", "𐀉": "du",
    "𐀊": "ja", "𐀋": "je",
    "𐀐": "ka", "𐀑": "ke", "𐀒": "ki", "𐀓": "ko", "𐀔": "ku",
    "𐀕": "ma", "𐀖": "me", "𐀗": "mi", "𐀘": "mo", "𐀙": "mu",
    "𐀚": "na", "𐀛": "ni", "𐀜": "no", "𐀝": "nu",
    "𐀞": "pa", "𐀟": "pe", "𐀠": "pi", "𐀡": "po", "𐀢": "pu",
    "𐀣": "qa", "𐀤": "qe", "𐀥": "qi", "𐀦": "qo",
    "𐀨": "ra", "𐀩": "re", "𐀪": "ri", "𐀫": "ro", "𐀬": "ru",
    "𐀭": "sa", "𐀮": "se", "𐀯": "si", "𐀰": "so", "𐀱": "su",
    "𐀲": "ta", "𐀳": "te", "𐀴": "ti", "𐀵": "to", "𐀶": "tu",
    "𐀷": "wa", "𐀸": "we", "𐀹": "wi", "𐀺": "wo",
    "𐀼": "za", "𐀽": "ze", "𐀾": "zo", "𐀿": "zu",
    # Variant signs
    "𐁀": "ha", "𐁁": "ai", "𐁂": "au", "𐁃": "dwe",
    "𐁄": "dwo", "𐁅": "nwa", "𐁆": "pu2", "𐁇": "pte",
    "𐁈": "ra3", "𐁉": "ro2", "𐁊": "ta2",
}

score_lines: list[str] = [
    "=== Learned Greek -> Linear B correspondences ===",
    "",
    "Best Linear B match for each common Greek character.",
    "LB_value shows the known phonetic reading of the Linear B sign.",
    "Only correspondences with score > 2.0 are shown (weaker ones are noise).",
    "",
    f"  {'Greek':>6}  {'LB sign':>8}  {'LB value':>9}  {'score':>8}",
    f"  {'-----':>6}  {'-------':>8}  {'--------':>9}  {'-----':>8}",
]

for g_char in sorted(best_per_greek, key=lambda c: -best_per_greek[c][1]):
    lb_char, score = best_per_greek[g_char]
    if score < 2.0:
        continue
    lb_val = LB_VALUES.get(lb_char, "?")
    score_lines.append(f"  {g_char:>6}  {lb_char:>8}  {lb_val:>9}  {score:>+8.3f}")

top_file = OUTPUT_DIR / "neurodecipher_top_scores.txt"
top_file.write_text("\n".join(score_lines), encoding="utf-8")
print(f"Learned correspondences written to {top_file}")
