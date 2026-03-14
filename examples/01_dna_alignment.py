"""Basic DNA sequence alignment with malign.

Demonstrates: align(), tabulate_alms(), DNA_MATRIX, k-best alignments,
method comparison (dumb / anw / yenksp).

No external data files needed -- sequences are hardcoded.
"""

from pathlib import Path

from malign import align, tabulate_alms
from malign.utils import DNA_MATRIX

OUTPUT_DIR = Path(__file__).parent / "output"
OUTPUT_DIR.mkdir(exist_ok=True)

# --- sequences ---
seqs = [
    ("GATTACA", "ATAC"),
    ("AGTACGCA", "TATGC"),
    ("GCATGCG", "GATTACA", "GCATGCG"),  # 3-way
]

lines: list[str] = []

# 1. Default identity matrix
lines.append("=== Pairwise alignment with default identity matrix ===")
alms = align(list(seqs[0]), k=1)
lines.append(tabulate_alms(alms))
lines.append("")

# 2. DNA scoring matrix
lines.append("=== Pairwise alignment with DNA_MATRIX ===")
alms = align(list(seqs[0]), matrix=DNA_MATRIX, k=1)
lines.append(tabulate_alms(alms))
lines.append("")

# 3. k-best alignments
lines.append("=== Top-4 alignments (yenksp, DNA_MATRIX) ===")
alms = align(list(seqs[0]), method="yenksp", matrix=DNA_MATRIX, k=4)
lines.append(tabulate_alms(alms))
lines.append("")

# 4. Method comparison
lines.append("=== Method comparison on second pair ===")
for method in ("dumb", "anw", "yenksp"):
    alms = align(list(seqs[1]), method=method, matrix=DNA_MATRIX, k=1)
    lines.append(f"--- {method} ---")
    lines.append(tabulate_alms(alms))
    lines.append("")

# 5. Three-way alignment
lines.append("=== Three-way alignment ===")
alms = align(list(seqs[2]), k=1)
lines.append(tabulate_alms(alms))

output = "\n".join(lines)
outfile = OUTPUT_DIR / "dna_alignments.txt"
outfile.write_text(output, encoding="utf-8")
print(output)
print(f"\nWritten to {outfile}")
