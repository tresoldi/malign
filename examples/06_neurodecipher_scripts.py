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
import math
from collections import Counter
from pathlib import Path

from malign import align, bootstrap_matrix

DATA_DIR = Path(__file__).parent / "data"
OUTPUT_DIR = Path(__file__).parent / "output"
OUTPUT_DIR.mkdir(exist_ok=True)

# Known Linear B phonetic values (syllabary readings)
LB_VALUES = {
    "\U00010000": "a", "\U00010001": "e", "\U00010002": "i",
    "\U00010003": "o", "\U00010004": "u",
    "\U00010005": "da", "\U00010006": "de", "\U00010007": "di",
    "\U00010008": "do", "\U00010009": "du",
    "\U0001000A": "ja", "\U0001000B": "je",
    "\U0001000D": "jo", "\U0001000F": "ka2",
    "\U00010010": "ka", "\U00010011": "ke", "\U00010012": "ki",
    "\U00010013": "ko", "\U00010014": "ku",
    "\U00010015": "ma", "\U00010016": "me", "\U00010017": "mi",
    "\U00010018": "mo", "\U00010019": "mu",
    "\U0001001A": "na", "\U0001001B": "ni", "\U0001001C": "no",
    "\U0001001D": "nu",
    "\U0001001E": "pa", "\U0001001F": "pe", "\U00010020": "pi",
    "\U00010021": "po", "\U00010022": "pu",
    "\U00010023": "qa", "\U00010024": "qe", "\U00010025": "qi",
    "\U00010026": "qo",
    "\U00010028": "ra", "\U00010029": "re", "\U0001002A": "ri",
    "\U0001002B": "ro", "\U0001002C": "ru",
    "\U0001002D": "sa", "\U0001002E": "se", "\U0001002F": "si",
    "\U00010030": "so", "\U00010031": "su",
    "\U00010032": "ta", "\U00010033": "te", "\U00010034": "ti",
    "\U00010035": "to", "\U00010036": "tu",
    "\U00010037": "wa", "\U00010038": "we", "\U00010039": "wi",
    "\U0001003A": "wo",
    "\U0001003C": "za", "\U0001003D": "ze", "\U0001003E": "zo",
    "\U0001003F": "zu",
    "\U00010040": "a2", "\U00010041": "a3", "\U00010042": "au",
    "\U00010043": "dwe", "\U00010044": "dwo", "\U00010045": "nwa",
    "\U00010046": "pu2", "\U00010047": "pte", "\U00010048": "ra2",
    "\U00010049": "ra3", "\U0001004A": "ro2", "\U0001004B": "ta2",
}


def lb_reading(sign: str) -> str:
    """Return the phonetic reading of a Linear B sign, or the sign itself."""
    if sign == "-":
        return "-"
    return LB_VALUES.get(sign, sign)


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
    alignment_method="yenksp",
    verbose=True,
)

# --- Step 2: Annotated alignments ---
# Select diverse, well-proportioned pairs (Greek/LB length ratio near 2.0)
scored_pairs: list[tuple[float, list[str], list[str]]] = []
for greek, linearb in all_pairs:
    if len(linearb) < 2:
        continue
    ratio = len(greek) / len(linearb)
    score = abs(ratio - 2.0)  # prefer ratio close to 2
    scored_pairs.append((score, greek, linearb))
scored_pairs.sort()

# Deduplicate by Greek word, pick diverse initial letters
seen_words: set[str] = set()
selected: list[tuple[list[str], list[str]]] = []
for _, greek, linearb in scored_pairs:
    word = "".join(greek)
    if word in seen_words:
        continue
    seen_words.add(word)
    selected.append((greek, linearb))
    if len(selected) >= 30:
        break

lines: list[str] = [
    "=== Greek (alphabet) -> Linear B (syllabary) ===",
    "",
    "Each alignment shows three rows:",
    "  Greek:   the alphabetic Greek characters",
    "  Linear B: the syllabary signs",
    "  Reading: phonetic reading of each Linear B sign",
    "",
    "Since each Linear B sign encodes a CV syllable, it 'absorbs'",
    "a Greek consonant+vowel pair, producing structural gaps.",
    "",
    "--- Identity matrix (no learning) ---",
    "",
]

for greek, linearb in selected[:10]:
    alms = align([greek, linearb], k=1)
    if not alms:
        continue
    a = alms[0]
    reading = [lb_reading(str(s)) for s in a.seqs[1]]
    lines.append(f"  {''.join(greek)} -> {''.join(linearb)}")
    lines.append(f"  Greek:    {' '.join(str(s) for s in a.seqs[0])}")
    lines.append(f"  Linear B: {' '.join(str(s) for s in a.seqs[1])}")
    lines.append(f"  Reading:  {' '.join(reading)}")
    lines.append(f"  score: {a.score:.2f}")
    lines.append("")

lines.append("--- Learned matrix (after bootstrap) ---")
lines.append("")

for greek, linearb in selected[:10]:
    alms = align([greek, linearb], matrix=learned, k=1)
    if not alms:
        continue
    a = alms[0]
    reading = [lb_reading(str(s)) for s in a.seqs[1]]
    lines.append(f"  {''.join(greek)} -> {''.join(linearb)}")
    lines.append(f"  Greek:    {' '.join(str(s) for s in a.seqs[0])}")
    lines.append(f"  Linear B: {' '.join(str(s) for s in a.seqs[1])}")
    lines.append(f"  Reading:  {' '.join(reading)}")
    lines.append(f"  score: {a.score:.2f}")
    lines.append("")

# Additional learned alignments (diverse set)
lines.append("--- More learned alignments ---")
lines.append("")

for greek, linearb in selected[10:30]:
    alms = align([greek, linearb], matrix=learned, k=1)
    if not alms:
        continue
    a = alms[0]
    reading = [lb_reading(str(s)) for s in a.seqs[1]]
    lines.append(f"  {''.join(greek)} -> {''.join(linearb)}")
    lines.append(f"  Greek:    {' '.join(str(s) for s in a.seqs[0])}")
    lines.append(f"  Linear B: {' '.join(str(s) for s in a.seqs[1])}")
    lines.append(f"  Reading:  {' '.join(reading)}")
    lines.append(f"  score: {a.score:.2f}")
    lines.append("")

alm_file = OUTPUT_DIR / "neurodecipher_alignments.txt"
alm_file.write_text("\n".join(lines), encoding="utf-8")
print(f"Alignments written to {alm_file}")

# --- Step 3: Learned correspondences ---
# Count co-occurrences: how often does each Greek char appear in a word
# alongside each Linear B sign?
cooccur: Counter[tuple[str, str]] = Counter()
greek_total: Counter[str] = Counter()
for greek, linearb in all_pairs:
    g_set = set(greek)
    lb_set = set(linearb)
    for g in g_set:
        greek_total[g] += 1
        for lb in lb_set:
            cooccur[(g, lb)] += 1

# For each Greek character, find the best Linear B matches from the learned
# matrix, but only among pairs that actually co-occur in the data.
gap = learned.gap


def best_matches(g_char: str, top_n: int = 3) -> list[tuple[str, float, int]]:
    """Return top-N (lb_sign, score, count) for a Greek char, filtered to co-occur > 0.

    Ranks by score * log(count + 1) to balance learned affinity with
    co-occurrence evidence -- pure score picks rare noise, pure count
    ignores what the matrix learned.
    """
    candidates: list[tuple[float, float, str, int]] = []
    for pair, score in learned.scores.items():
        if pair[0] == g_char and gap not in pair:
            count = cooccur.get((g_char, pair[1]), 0)
            if count > 0 and pair[1] in LB_VALUES:
                rank = score * math.log1p(count)
                candidates.append((rank, score, pair[1], count))
    candidates.sort(reverse=True)
    return [(lb, sc, ct) for _, sc, lb, ct in candidates[:top_n]]


score_lines: list[str] = [
    "=== Learned Greek -> Linear B correspondences ===",
    "",
    "For each Greek character: top Linear B signs from the learned matrix,",
    "filtered to pairs that actually co-occur in the data.",
    "(Phonetically correct matches are marked with *)",
    "",
    f"  {'Greek':>6}  {'LB sign':<6}  {'reading':<8}  {'score':>8}  {'co-occur':>8}",
    f"  {'-----':>6}  {'------':<6}  {'-------':<8}  {'-----':>8}  {'--------':>8}",
]

# Known correct Greek -> Linear B consonant correspondences
# (Greek consonant maps to LB CV-sign starting with same consonant)
CORRECT = {
    "δ": {"da", "de", "di", "do", "du"},
    "ζ": {"za", "ze", "zo", "zu"},
    "θ": {"ta", "ta2", "te", "ti", "to", "tu"},  # θ > t in Mycenaean
    "κ": {"ka", "ka2", "ke", "ki", "ko", "ku"},
    "λ": {"ra", "re", "ri", "ro", "ru", "ra2", "ra3", "ro2"},  # l/r merge in LB
    "μ": {"ma", "me", "mi", "mo", "mu"},
    "ν": {"na", "ni", "no", "nu", "nwa"},
    "π": {"pa", "pe", "pi", "po", "pu", "pu2", "pte"},
    "ρ": {"ra", "re", "ri", "ro", "ru", "ra2", "ra3", "ro2"},
    "σ": {"sa", "se", "si", "so", "su"},
    "ς": {"sa", "se", "si", "so", "su"},
    "τ": {"ta", "ta2", "te", "ti", "to", "tu"},
    "φ": {"pa", "pe", "pi", "po", "pu", "pu2"},  # φ > p in Mycenaean
    "χ": {"ka", "ka2", "ke", "ki", "ko", "ku"},  # χ > k in Mycenaean
    "α": {"a"},
    "ε": {"e"},
    "ι": {"i"},
    "ο": {"o"},
    "υ": {"u"},
}

greek_consonants = "δζθκλμνπρσςτφχ"  # only those with enough data

for g_char in greek_consonants:
    matches = best_matches(g_char)
    if not matches:
        continue
    for i, (lb, sc, ct) in enumerate(matches):
        reading = lb_reading(lb)
        correct_set = CORRECT.get(g_char, set())
        mark = " *" if reading in correct_set else ""
        prefix = g_char if i == 0 else ""
        score_lines.append(
            f"  {prefix:>6}  {lb:<6}  {reading:<8}  {sc:>+8.3f}  {ct:>8}{mark}"
        )
    if len(matches) > 0:
        score_lines.append("")

# Vowels
score_lines.append("  Vowels:")
for g_char in "αεηιου":
    matches = best_matches(g_char)
    if not matches:
        continue
    for i, (lb, sc, ct) in enumerate(matches):
        reading = lb_reading(lb)
        correct_set = CORRECT.get(g_char, set())
        mark = " *" if reading in correct_set else ""
        prefix = g_char if i == 0 else ""
        score_lines.append(
            f"  {prefix:>6}  {lb:<6}  {reading:<8}  {sc:>+8.3f}  {ct:>8}{mark}"
        )
    if len(matches) > 0:
        score_lines.append("")

# Summary
n_correct = 0
n_total = 0
for g_char in greek_consonants + "αεηιου":
    matches = best_matches(g_char, top_n=1)
    if matches:
        n_total += 1
        lb, _, _ = matches[0]
        reading = lb_reading(lb)
        if reading in CORRECT.get(g_char, set()):
            n_correct += 1

score_lines.append(
    f"  Top-1 accuracy: {n_correct}/{n_total} correct"
    f" ({100 * n_correct / n_total:.0f}%)"
)

top_file = OUTPUT_DIR / "neurodecipher_top_scores.txt"
top_file.write_text("\n".join(score_lines), encoding="utf-8")
print(f"Learned correspondences written to {top_file}")
print(f"  Top-1 accuracy: {n_correct}/{n_total}")
