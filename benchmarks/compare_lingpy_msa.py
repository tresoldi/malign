"""MSA benchmark comparing malign against LingPy on gold-standard cognate sets.

Unlike the pairwise benchmark, this evaluates true multiple sequence alignment
(3+ sequences per cognate set). malign uses exact N-dimensional alignment for
small groups (≤4 seqs) and UPGMA-guided progressive for larger ones. LingPy
uses UPGMA-guided progressive alignment with SCA sound classes.

Gold alignments come from BDPA datasets and NorthEuralex (MSA-derived _REF
columns). All-gap columns are stripped before evaluation.

Datasets: all 12 BDPA files + NorthEuralex Italic (italic_by_cogid.tsv).

Usage:
  pip install malign[lingpy,features]
  python benchmarks/compare_lingpy_msa.py
"""

import csv
import itertools
import logging
import signal
import sys
import time
from pathlib import Path

from tabulate import tabulate

# Suppress LingPy's verbose INFO logging
logging.getLogger("lingpy").setLevel(logging.WARNING)

from malign import (
    Alignment,
    align,
    alignment_accuracy,
    alignment_f1,
    alignment_precision_recall,
    bootstrap_matrix,
    strip_common_gaps,
)
from malign.anw import nw_align
from malign.ndim_common import should_use_ndim
from malign.progressive import upgma_progressive_align
from malign.scoring_matrix import ScoringMatrix

try:
    from lingpy import Multiple
except ImportError:
    sys.exit(
        "lingpy not installed. Install with:\n"
        "  pip install malign[lingpy]\n"
        "or:\n"
        "  pip install lingpy"
    )

# --- Paths ---
DATA_DIR = Path(__file__).resolve().parent.parent / "examples" / "data"
OUTPUT_DIR = Path(__file__).resolve().parent / "output"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

GAP = "-"

# Max sequences per cognate set — malign progressive scales ~10x per seq beyond 5
MAX_SEQS = 8

# Per-alignment timeout in seconds
ALIGN_TIMEOUT = 5


class _AlignTimeout(Exception):
    pass


def _timeout_handler(signum, frame):
    raise _AlignTimeout

# All BDPA datasets: (label, filename)
BDPA_DATASETS = [
    ("BDPA Romance", "romance.tsv"),
    ("BDPA Slavic", "slavic.tsv"),
    ("BDPA Bai", "bai.tsv"),
    ("BDPA Germanic", "germanic.tsv"),
    ("BDPA Sinitic", "sinitic.tsv"),
    ("BDPA Andean", "andean.tsv"),
    ("BDPA French", "french.tsv"),
    ("BDPA Japanese", "japanese.tsv"),
    ("BDPA Norwegian", "norwegian.tsv"),
    ("BDPA Ob-Ugrian", "ob-ugrian.tsv"),
    ("BDPA Bulgarian", "bulgarian.tsv"),
    ("BDPA Dutch", "dutch.tsv"),
]

# Type for a cognate set
CognateSet = dict[str, object]


# ============================================================================
# Data Loaders
# ============================================================================


def load_bdpa_msa(filename: str, min_seqs: int = 3, max_seqs: int = MAX_SEQS) -> list[CognateSet]:
    """Load a BDPA TSV file as MSA cognate sets."""
    rows: list[dict[str, str]] = []
    with open(DATA_DIR / "bdpa" / filename, encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            rows.append(row)

    header_langs = []
    for col in rows[0]:
        if col.endswith("_SEQ"):
            header_langs.append(col.removesuffix("_SEQ"))

    cognate_sets: list[CognateSet] = []
    for row in rows:
        gloss = row["ID"]
        languages = []
        sequences = []
        references = []

        for lang in header_langs:
            ref = row.get(f"{lang}_REF", "").strip()
            seq = row.get(f"{lang}_SEQ", "").strip()
            if ref and seq:
                languages.append(lang)
                sequences.append(seq.split())
                references.append(ref.split())

        if len(languages) < min_seqs:
            continue

        # Cap at max_seqs languages
        if len(languages) > max_seqs:
            languages = languages[:max_seqs]
            sequences = sequences[:max_seqs]
            references = references[:max_seqs]

        # Validate: all references must have same length (valid MSA)
        ref_lens = {len(r) for r in references}
        if len(ref_lens) != 1:
            continue

        # Strip all-gap columns from gold MSA
        gold = strip_common_gaps(Alignment(seqs=references, score=0.0))

        cognate_sets.append({
            "id": gloss,
            "languages": languages,
            "sequences": sequences,
            "references": [list(s) for s in gold.seqs],
            "num_seqs": len(languages),
        })

    return cognate_sets


def load_northeuralex_msa(min_seqs: int = 3, max_seqs: int = MAX_SEQS) -> list[CognateSet]:
    """Load NorthEuralex italic_by_cogid.tsv as MSA cognate sets."""
    rows: list[dict[str, str]] = []
    with open(DATA_DIR / "northeuralex" / "italic_by_cogid.tsv", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            rows.append(row)

    header_langs = []
    for col in rows[0]:
        if col.endswith("_SEQ"):
            header_langs.append(col.removesuffix("_SEQ"))

    cognate_sets: list[CognateSet] = []
    for row in rows:
        cogid = row.get("COGID", "")
        languages = []
        sequences = []
        references = []

        for lang in header_langs:
            ref = row.get(f"{lang}_REF", "").strip()
            seq = row.get(f"{lang}_SEQ", "").strip()
            if ref and seq:
                languages.append(lang)
                sequences.append([s for s in seq.split() if s != "+"])
                references.append([s for s in ref.split() if s != "+"])

        if len(languages) < min_seqs:
            continue

        # Cap at max_seqs languages
        if len(languages) > max_seqs:
            languages = languages[:max_seqs]
            sequences = sequences[:max_seqs]
            references = references[:max_seqs]

        # Validate: all references must have same length (valid MSA)
        ref_lens = {len(r) for r in references}
        if len(ref_lens) != 1:
            continue

        # Strip all-gap columns from gold MSA
        gold = strip_common_gaps(Alignment(seqs=references, score=0.0))

        cognate_sets.append({
            "id": cogid,
            "languages": languages,
            "sequences": sequences,
            "references": [list(s) for s in gold.seqs],
            "num_seqs": len(languages),
        })

    return cognate_sets


# ============================================================================
# Alignment Wrappers
# ============================================================================


def align_lingpy_msa(sequences: list[list[str]]) -> Alignment | None:
    """Align multiple sequences using LingPy's progressive alignment."""
    try:
        strs = [" ".join(seq) for seq in sequences]
        m = Multiple(strs)
        m.prog_align(model="sca", mode="global")
        return Alignment(seqs=[list(row) for row in m.alm_matrix], score=0.0)
    except Exception:
        return None


def align_malign_auto(sequences: list[list[str]], matrix=None) -> Alignment | None:
    """Align using malign's auto-dispatch (ndim for small, progressive for large)."""
    try:
        alms = align(sequences, k=1, matrix=matrix)
        return alms[0] if alms else None
    except Exception:
        return None


def align_malign_progressive(sequences: list[list[str]], matrix=None) -> Alignment | None:
    """Force progressive alignment (bypass ndim dispatch)."""
    try:
        seqs = [list(s) for s in sequences]
        if matrix is None:
            matrix = ScoringMatrix.from_sequences(seqs, match=+1, gap_score=-1)
        alms = upgma_progressive_align(seqs, matrix, pw_func=nw_align, k=1)
        return alms[0] if alms else None
    except Exception:
        return None


# ============================================================================
# Evaluation
# ============================================================================


def classify_dispatch(cog: CognateSet) -> str:
    """Classify whether a cognate set would use ndim or progressive."""
    n = cog["num_seqs"]
    seq_lengths = [len(s) for s in cog["sequences"]]
    return "ndim" if should_use_ndim(n, seq_lengths, "anw") else "progressive"


def evaluate_msa(
    cognate_sets: list[CognateSet],
    align_fn,
    method_label: str,
    dataset_label: str,
) -> tuple[list[dict], float]:
    """Evaluate MSA alignment function on cognate sets."""
    results = []
    t0 = time.perf_counter()

    for idx, cog in enumerate(cognate_sets):
        sequences = cog["sequences"]
        references = cog["references"]
        dispatch = classify_dispatch(cog)

        # Progress indicator
        if (idx + 1) % 10 == 0 or idx == len(cognate_sets) - 1:
            print(f"    {method_label}: {idx + 1}/{len(cognate_sets)}", flush=True)

        try:
            old_handler = signal.signal(signal.SIGALRM, _timeout_handler)
            signal.alarm(ALIGN_TIMEOUT)
            predicted = align_fn(sequences)
            signal.alarm(0)
            signal.signal(signal.SIGALRM, old_handler)
        except _AlignTimeout:
            signal.alarm(0)
            signal.signal(signal.SIGALRM, old_handler)
            predicted = None

        if predicted is None:
            results.append({
                "Dataset": dataset_label,
                "ID": cog["id"],
                "NumSeqs": cog["num_seqs"],
                "Dispatch": dispatch,
                "Method": method_label,
                "ColScore": None,
                "Precision": None,
                "Recall": None,
                "F1": None,
                "Status": "error",
            })
            continue

        # Strip all-gap columns from predicted
        predicted = strip_common_gaps(predicted)
        gold = Alignment(seqs=references, score=0.0)

        # Verify sequence count matches
        if len(predicted.seqs) != len(gold.seqs):
            results.append({
                "Dataset": dataset_label,
                "ID": cog["id"],
                "NumSeqs": cog["num_seqs"],
                "Dispatch": dispatch,
                "Method": method_label,
                "ColScore": None,
                "Precision": None,
                "Recall": None,
                "F1": None,
                "Status": "error",
            })
            continue

        # Column score: requires same alignment length
        col_score = None
        if len(predicted.seqs[0]) == len(gold.seqs[0]):
            col_score = alignment_accuracy(predicted, gold)

        # Pair-based metrics: work across different lengths
        try:
            prec, rec = alignment_precision_recall(predicted, gold)
            f1 = alignment_f1(predicted, gold)
        except Exception:
            prec, rec, f1 = None, None, None

        status = "ok" if col_score is not None else "len_mismatch"

        results.append({
            "Dataset": dataset_label,
            "ID": cog["id"],
            "NumSeqs": cog["num_seqs"],
            "Dispatch": dispatch,
            "Method": method_label,
            "ColScore": col_score,
            "Precision": prec,
            "Recall": rec,
            "F1": f1,
            "Status": status,
        })

    elapsed = time.perf_counter() - t0
    return results, elapsed


def summarize_msa(results: list[dict], elapsed: float) -> dict:
    """Compute aggregate statistics from MSA result dicts."""
    has_col = [r for r in results if r["ColScore"] is not None]
    has_pr = [r for r in results if r["Precision"] is not None]
    n_err = sum(1 for r in results if r["Status"] == "error")
    n_mismatch = sum(1 for r in results if r["Status"] == "len_mismatch")

    return {
        "Avg_ColScore": sum(r["ColScore"] for r in has_col) / len(has_col) if has_col else 0.0,
        "Avg_Precision": sum(r["Precision"] for r in has_pr) / len(has_pr) if has_pr else 0.0,
        "Avg_Recall": sum(r["Recall"] for r in has_pr) / len(has_pr) if has_pr else 0.0,
        "Avg_F1": sum(r["F1"] for r in has_pr) / len(has_pr) if has_pr else 0.0,
        "Num_Sets": len(has_col),
        "Num_LenMismatch": n_mismatch,
        "Num_Errors": n_err,
        "Avg_NumSeqs": (
            sum(r["NumSeqs"] for r in results) / len(results) if results else 0
        ),
        "Runtime_Seconds": elapsed,
    }


# ============================================================================
# Helpers
# ============================================================================


def collect_all_symbols(cognate_sets: list[CognateSet]) -> set[str]:
    """Collect all unique symbols across all sequences in cognate sets."""
    symbols: set[str] = set()
    for cog in cognate_sets:
        for seq in cog["sequences"]:
            symbols.update(seq)
    return symbols


def extract_pairs(cognate_sets: list[CognateSet], max_pairs: int = 200) -> list[tuple]:
    """Extract pairwise sequence pairs from cognate sets for bootstrap."""
    pairs = []
    for cog in cognate_sets:
        seqs = cog["sequences"]
        for i, j in itertools.combinations(range(len(seqs)), 2):
            pairs.append((seqs[i], seqs[j]))
    # Sample if too many
    if len(pairs) > max_pairs:
        import random
        rng = random.Random(42)
        pairs = rng.sample(pairs, max_pairs)
    return pairs


def run_dataset(
    dataset_label: str,
    cognate_sets: list[CognateSet],
    methods: list[tuple[str, object]],
) -> tuple[list[dict], list[dict]]:
    """Run all methods on a dataset and return (all_details, summary_rows)."""
    all_details = []
    summary_rows = []
    table_rows = []

    for method_label, align_fn in methods:
        results, elapsed = evaluate_msa(cognate_sets, align_fn, method_label, dataset_label)
        all_details.extend(results)

        stats = summarize_msa(results, elapsed)
        summary_rows.append({
            "Dataset": dataset_label,
            "Method": method_label,
            **stats,
        })

        table_rows.append([
            method_label,
            f"{stats['Avg_ColScore']:.3f}",
            f"{stats['Avg_F1']:.3f}",
            stats["Num_Sets"],
            stats["Num_LenMismatch"],
            stats["Num_Errors"],
            f"{stats['Avg_NumSeqs']:.1f}",
            f"{stats['Runtime_Seconds']:.1f}s",
        ])

    headers = ["Method", "Col", "F1", "Sets", "LenMM", "Err", "AvgN", "Time"]
    print(f"\n{dataset_label}")
    print(tabulate(table_rows, headers=headers, tablefmt="simple"))

    return all_details, summary_rows


# ============================================================================
# Main
# ============================================================================


def main():
    print("=== malign vs LingPy MSA Benchmark ===")
    print(f"    MAX_SEQS={MAX_SEQS}, TIMEOUT={ALIGN_TIMEOUT}s")

    all_details: list[dict] = []
    all_summaries: list[dict] = []

    # Collect all datasets
    datasets: list[tuple[str, list[CognateSet]]] = []

    for label, filename in BDPA_DATASETS:
        cogs = load_bdpa_msa(filename, min_seqs=3)
        if cogs:
            datasets.append((label, cogs))

    nel_cogs = load_northeuralex_msa(min_seqs=3)
    if nel_cogs:
        datasets.append(("NorthEuralex Italic", nel_cogs))

    for dataset_label, cognate_sets in datasets:
        n_sets = len(cognate_sets)
        sizes = [c["num_seqs"] for c in cognate_sets]
        n_ndim = sum(1 for c in cognate_sets if classify_dispatch(c) == "ndim")
        print(f"\n{'=' * 60}")
        print(f"{dataset_label}: {n_sets} cognate sets, "
              f"{min(sizes)}-{max(sizes)} seqs/set (avg {sum(sizes)/len(sizes):.1f}), "
              f"{n_ndim} ndim-eligible")

        # Build SCA matrix
        symbols = collect_all_symbols(cognate_sets)
        sca_matrix = ScoringMatrix.from_lingpy([sorted(symbols)], model="sca")

        # Bootstrap from pairwise extractions
        pairs = extract_pairs(cognate_sets, max_pairs=200)
        if pairs:
            print(f"  Bootstrapping from {len(pairs)} pairs...")
            boot_matrix = bootstrap_matrix(pairs, max_iter=5, verbose=False)
            boot_sca_matrix = bootstrap_matrix(
                pairs, max_iter=5, prior_matrix=sca_matrix, prior_weight=0.3, verbose=False,
            )
            print("  Bootstrap done.")
        else:
            boot_matrix = None
            boot_sca_matrix = None

        # --- Methods ---
        # Note: identity matrix is omitted — with N-domain matrices the progressive
        # alignment attempts ndim subtrees, making it orders of magnitude slower
        # than pairwise-matrix methods with no benefit for quality.
        methods = [
            ("LingPy prog", lambda seqs: align_lingpy_msa(seqs)),
            ("malign auto (SCA)", lambda seqs, m=sca_matrix: align_malign_auto(seqs, matrix=m)),
        ]

        if boot_matrix is not None:
            methods.append(
                ("malign auto (boot)",
                 lambda seqs, m=boot_matrix: align_malign_auto(seqs, matrix=m))
            )
            methods.append(
                ("malign auto (SCA+boot)",
                 lambda seqs, m=boot_sca_matrix: align_malign_auto(seqs, matrix=m))
            )

        # For ndim-eligible sets, also run forced progressive to compare dispatch strategies
        if n_ndim > 0:
            methods.append(
                ("malign prog (SCA)",
                 lambda seqs, m=sca_matrix: align_malign_progressive(seqs, matrix=m))
            )

        details, summaries = run_dataset(dataset_label, cognate_sets, methods)
        all_details.extend(details)
        all_summaries.extend(summaries)

    # --- Write output files ---
    summary_file = OUTPUT_DIR / "compare_lingpy_msa_summary.tsv"
    summary_fields = [
        "Dataset", "Method", "Avg_ColScore", "Avg_Precision", "Avg_Recall",
        "Avg_F1", "Num_Sets", "Num_LenMismatch", "Num_Errors",
        "Avg_NumSeqs", "Runtime_Seconds",
    ]
    with open(summary_file, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=summary_fields, delimiter="\t")
        writer.writeheader()
        for row in all_summaries:
            writer.writerow({k: f"{v:.4f}" if isinstance(v, float) else v for k, v in row.items()})
    print(f"\nSummary written to {summary_file}")

    detail_file = OUTPUT_DIR / "compare_lingpy_msa_detail.tsv"
    detail_fields = [
        "Dataset", "ID", "NumSeqs", "Dispatch", "Method",
        "ColScore", "Precision", "Recall", "F1", "Status",
    ]
    with open(detail_file, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=detail_fields, delimiter="\t")
        writer.writeheader()
        for row in all_details:
            out = {}
            for k in detail_fields:
                v = row[k]
                out[k] = f"{v:.4f}" if isinstance(v, float) else (v if v is not None else "")
            writer.writerow(out)
    print(f"Details written to {detail_file}")


if __name__ == "__main__":
    main()
