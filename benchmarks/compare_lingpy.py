"""Benchmark comparing malign alignment quality against LingPy on gold-standard datasets.

LingPy uses static scoring matrices (SCA/DOLGO/ASJP) while malign can learn
matrices from data via bootstrap/EM. Evaluation uses malign's own metrics for
apples-to-apples comparison.

Datasets:
  - BDPA Romance (8 languages, ~840 pairwise comparisons)
  - NorthEuralex Italic (500 sampled pairs)
  - Lexibank savelyevturkic (Turkish-Azeri, ~200 pairs)

Usage:
  pip install malign[lingpy]
  python benchmarks/compare_lingpy.py
"""

import csv
import itertools
import logging
import random
import sys
import time
from collections import defaultdict
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
    learn_matrix,
)

try:
    from lingpy import Pairwise
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

# Type alias for pair data
PairData = dict[str, object]


# ============================================================================
# B. Data Loaders
# ============================================================================


def load_bdpa() -> list[PairData]:
    """Load BDPA Romance cognate pairs with gold alignments."""
    rows: list[dict[str, str]] = []
    with open(DATA_DIR / "bdpa" / "romance.tsv", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            rows.append(row)

    # Extract language names from _SEQ column suffixes
    header_langs = []
    for col in rows[0]:
        if col.endswith("_SEQ"):
            header_langs.append(col.removesuffix("_SEQ"))

    pairs: list[PairData] = []
    for row in rows:
        gloss = row["ID"]

        # Collect languages with data in this row
        available = []
        for lang in header_langs:
            ref = row.get(f"{lang}_REF", "").strip()
            seq = row.get(f"{lang}_SEQ", "").strip()
            if ref and seq:
                available.append(lang)

        # Generate all pairwise combinations
        for lang_a, lang_b in itertools.combinations(available, 2):
            pairs.append({
                "id": gloss,
                "lang_a": lang_a,
                "lang_b": lang_b,
                "seq_a": row[f"{lang_a}_SEQ"].split(),
                "seq_b": row[f"{lang_b}_SEQ"].split(),
                "ref_a": row[f"{lang_a}_REF"].split(),
                "ref_b": row[f"{lang_b}_REF"].split(),
            })

    return pairs


def load_northeuralex(n_sample: int = 500, seed: int = 42) -> list[PairData]:
    """Load NorthEuralex Italic cognate pairs, sampled for evaluation."""
    rows: list[dict[str, str]] = []
    with open(DATA_DIR / "northeuralex" / "italic_by_cogid.tsv", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            rows.append(row)

    # Extract language names from _SEQ column suffixes
    header_langs = []
    for col in rows[0]:
        if col.endswith("_SEQ"):
            header_langs.append(col.removesuffix("_SEQ"))

    pairs: list[PairData] = []
    for row in rows:
        cogid = row.get("COGID", "")

        available = []
        for lang in header_langs:
            ref = row.get(f"{lang}_REF", "").strip()
            seq = row.get(f"{lang}_SEQ", "").strip()
            if ref and seq:
                available.append(lang)

        for lang_a, lang_b in itertools.combinations(available, 2):
            # Filter morpheme boundary markers
            seq_a = [s for s in row[f"{lang_a}_SEQ"].split() if s != "+"]
            seq_b = [s for s in row[f"{lang_b}_SEQ"].split() if s != "+"]
            ref_a = [s for s in row[f"{lang_a}_REF"].split() if s != "+"]
            ref_b = [s for s in row[f"{lang_b}_REF"].split() if s != "+"]

            if seq_a and seq_b and ref_a and ref_b:
                pairs.append({
                    "id": cogid,
                    "lang_a": lang_a,
                    "lang_b": lang_b,
                    "seq_a": seq_a,
                    "seq_b": seq_b,
                    "ref_a": ref_a,
                    "ref_b": ref_b,
                })

    # Sample for evaluation
    rng = random.Random(seed)
    if len(pairs) > n_sample:
        pairs = rng.sample(pairs, n_sample)

    return pairs


def load_lexibank() -> list[PairData]:
    """Load Lexibank savelyevturkic Turkish-Azeri cognate pairs."""
    dataset = "savelyevturkic"
    lang_a = "savelyevturkic_Turkish"
    lang_b = "savelyevturkic_Azeri"

    cognates: dict[str, dict[str, tuple[list[str], list[str]]]] = defaultdict(dict)

    with open(DATA_DIR / "lexibank" / "forms.tsv", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if row["Dataset"] != dataset:
                continue
            lang = row["Language_ID"]
            if lang not in (lang_a, lang_b):
                continue
            cog = row["Cognacy"].strip()
            seg = row["Segments"].strip()
            alm = row["Alignment"].strip()
            if cog and seg:
                cognates[cog][lang] = (seg.split(), alm.split() if alm else [])

    pairs: list[PairData] = []
    for cog_id, langs in cognates.items():
        if lang_a in langs and lang_b in langs:
            seg_a, alm_a = langs[lang_a]
            seg_b, alm_b = langs[lang_b]
            if alm_a and alm_b:
                pairs.append({
                    "id": cog_id,
                    "lang_a": "Turkish",
                    "lang_b": "Azeri",
                    "seq_a": seg_a,
                    "seq_b": seg_b,
                    "ref_a": alm_a,
                    "ref_b": alm_b,
                })

    return pairs


# ============================================================================
# C. Alignment Wrappers
# ============================================================================


def align_lingpy(seq_a: list[str], seq_b: list[str]) -> Alignment | None:
    """Align two sequences using LingPy's SCA model."""
    try:
        p = Pairwise(" ".join(seq_a), " ".join(seq_b))
        p.align(model="sca", mode="global")
        aligned_a, aligned_b, score = p.alignments[0]
        return Alignment(seqs=[list(aligned_a), list(aligned_b)], score=float(score))
    except Exception:
        return None


def align_malign(seq_a: list[str], seq_b: list[str], matrix=None) -> Alignment | None:
    """Align two sequences using malign."""
    alms = align([seq_a, seq_b], k=1, matrix=matrix)
    return alms[0] if alms else None


# ============================================================================
# D. Evaluation
# ============================================================================


def evaluate_pairs(
    pairs: list[PairData],
    align_fn,
    method_label: str,
    dataset_label: str,
) -> list[dict]:
    """Evaluate alignment function on pairs, returning per-pair results."""
    results = []

    t0 = time.perf_counter()
    for pair in pairs:
        seq_a = pair["seq_a"]
        seq_b = pair["seq_b"]
        ref_a = pair["ref_a"]
        ref_b = pair["ref_b"]

        predicted = align_fn(seq_a, seq_b)

        if predicted is None:
            results.append({
                "Dataset": dataset_label,
                "ID": pair["id"],
                "Lang_A": pair["lang_a"],
                "Lang_B": pair["lang_b"],
                "Method": method_label,
                "Accuracy": None,
                "Precision": None,
                "Recall": None,
                "F1": None,
                "Status": "error",
            })
            continue

        # Build gold alignment — skip if ref lengths don't match each other
        if len(ref_a) != len(ref_b):
            results.append({
                "Dataset": dataset_label,
                "ID": pair["id"],
                "Lang_A": pair["lang_a"],
                "Lang_B": pair["lang_b"],
                "Method": method_label,
                "Accuracy": None,
                "Precision": None,
                "Recall": None,
                "F1": None,
                "Status": "skipped",
            })
            continue

        gold = Alignment(seqs=[ref_a, ref_b], score=0.0)

        # Check length compatibility
        if len(predicted.seqs[0]) != len(gold.seqs[0]):
            results.append({
                "Dataset": dataset_label,
                "ID": pair["id"],
                "Lang_A": pair["lang_a"],
                "Lang_B": pair["lang_b"],
                "Method": method_label,
                "Accuracy": None,
                "Precision": None,
                "Recall": None,
                "F1": None,
                "Status": "skipped",
            })
            continue

        acc = alignment_accuracy(predicted, gold)
        prec, rec = alignment_precision_recall(predicted, gold)
        f1 = alignment_f1(predicted, gold)

        results.append({
            "Dataset": dataset_label,
            "ID": pair["id"],
            "Lang_A": pair["lang_a"],
            "Lang_B": pair["lang_b"],
            "Method": method_label,
            "Accuracy": acc,
            "Precision": prec,
            "Recall": rec,
            "F1": f1,
            "Status": "ok",
        })

    elapsed = time.perf_counter() - t0
    return results, elapsed


def summarize(results: list[dict], elapsed: float) -> dict:
    """Compute aggregate statistics from result dicts."""
    ok = [r for r in results if r["Status"] == "ok"]
    n_ok = len(ok)
    n_skip = sum(1 for r in results if r["Status"] == "skipped")
    n_err = sum(1 for r in results if r["Status"] == "error")

    if n_ok == 0:
        return {
            "Avg_Accuracy": 0.0,
            "Avg_Precision": 0.0,
            "Avg_Recall": 0.0,
            "Avg_F1": 0.0,
            "Num_Pairs": n_ok,
            "Num_Skipped": n_skip,
            "Num_Errors": n_err,
            "Runtime_Seconds": elapsed,
        }

    return {
        "Avg_Accuracy": sum(r["Accuracy"] for r in ok) / n_ok,
        "Avg_Precision": sum(r["Precision"] for r in ok) / n_ok,
        "Avg_Recall": sum(r["Recall"] for r in ok) / n_ok,
        "Avg_F1": sum(r["F1"] for r in ok) / n_ok,
        "Num_Pairs": n_ok,
        "Num_Skipped": n_skip,
        "Num_Errors": n_err,
        "Runtime_Seconds": elapsed,
    }


# ============================================================================
# E. Main
# ============================================================================


def run_dataset(
    dataset_label: str,
    pairs: list[PairData],
    methods: list[tuple[str, object]],
) -> tuple[list[dict], list[dict]]:
    """Run all methods on a dataset and return (all_details, summary_rows)."""
    all_details = []
    summary_rows = []
    table_rows = []

    for method_label, align_fn in methods:
        results, elapsed = evaluate_pairs(pairs, align_fn, method_label, dataset_label)
        all_details.extend(results)

        stats = summarize(results, elapsed)
        summary_rows.append({
            "Dataset": dataset_label,
            "Method": method_label,
            **stats,
        })

        table_rows.append([
            method_label,
            f"{stats['Avg_Accuracy']:.3f}",
            f"{stats['Avg_Precision']:.3f}",
            f"{stats['Avg_Recall']:.3f}",
            f"{stats['Avg_F1']:.3f}",
            stats["Num_Pairs"],
            stats["Num_Skipped"],
            f"{stats['Runtime_Seconds']:.1f}s",
        ])

    headers = ["Method", "Acc", "Prec", "Rec", "F1", "Pairs", "Skip", "Time"]
    print(f"\n{dataset_label}")
    print(tabulate(table_rows, headers=headers, tablefmt="simple"))

    return all_details, summary_rows


def main():
    print("=== malign vs LingPy Benchmark ===")

    all_details: list[dict] = []
    all_summaries: list[dict] = []

    # --- BDPA Romance ---
    print("\nLoading BDPA Romance data...")
    bdpa_pairs = load_bdpa()
    print(f"  {len(bdpa_pairs)} pairwise comparisons")

    print("  Bootstrapping matrix...")
    bdpa_bootstrap = bootstrap_matrix(
        [(p["seq_a"], p["seq_b"]) for p in bdpa_pairs],
        max_iter=15,
        verbose=False,
    )

    details, summaries = run_dataset(
        "BDPA Romance",
        bdpa_pairs,
        [
            ("LingPy SCA", lambda a, b: align_lingpy(a, b)),
            ("malign (identity)", lambda a, b: align_malign(a, b)),
            ("malign (bootstrap)", lambda a, b, m=bdpa_bootstrap: align_malign(a, b, matrix=m)),
        ],
    )
    all_details.extend(details)
    all_summaries.extend(summaries)

    # --- NorthEuralex Italic ---
    print("\nLoading NorthEuralex Italic data...")
    nel_pairs = load_northeuralex(n_sample=500, seed=42)
    print(f"  {len(nel_pairs)} pairwise comparisons (sampled)")

    print("  Bootstrapping matrix...")
    nel_bootstrap_pairs = [(p["seq_a"], p["seq_b"]) for p in nel_pairs[:200]]
    nel_bootstrap = bootstrap_matrix(
        nel_bootstrap_pairs,
        max_iter=15,
        verbose=False,
    )

    details, summaries = run_dataset(
        "NorthEuralex Italic",
        nel_pairs,
        [
            ("LingPy SCA", lambda a, b: align_lingpy(a, b)),
            ("malign (identity)", lambda a, b: align_malign(a, b)),
            ("malign (bootstrap)", lambda a, b, m=nel_bootstrap: align_malign(a, b, matrix=m)),
        ],
    )
    all_details.extend(details)
    all_summaries.extend(summaries)

    # --- Lexibank savelyevturkic ---
    print("\nLoading Lexibank savelyevturkic data...")
    lex_pairs = load_lexibank()
    print(f"  {len(lex_pairs)} pairwise comparisons")

    print("  Learning matrix...")
    lex_cognate_sets = [[p["seq_a"], p["seq_b"]] for p in lex_pairs]
    lex_learned = learn_matrix(
        lex_cognate_sets,
        max_iter=5,
        verbose=False,
    )

    details, summaries = run_dataset(
        "Lexibank Turkish-Azeri",
        lex_pairs,
        [
            ("LingPy SCA", lambda a, b: align_lingpy(a, b)),
            ("malign (identity)", lambda a, b: align_malign(a, b)),
            ("malign (learned)", lambda a, b, m=lex_learned: align_malign(a, b, matrix=m)),
        ],
    )
    all_details.extend(details)
    all_summaries.extend(summaries)

    # --- Write output files ---
    # Summary TSV
    summary_file = OUTPUT_DIR / "compare_lingpy_summary.tsv"
    summary_fields = [
        "Dataset", "Method", "Avg_Accuracy", "Avg_Precision", "Avg_Recall",
        "Avg_F1", "Num_Pairs", "Num_Skipped", "Num_Errors", "Runtime_Seconds",
    ]
    with open(summary_file, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=summary_fields, delimiter="\t")
        writer.writeheader()
        for row in all_summaries:
            writer.writerow({k: f"{v:.4f}" if isinstance(v, float) else v for k, v in row.items()})
    print(f"\nSummary written to {summary_file}")

    # Detail TSV
    detail_file = OUTPUT_DIR / "compare_lingpy_detail.tsv"
    detail_fields = [
        "Dataset", "ID", "Lang_A", "Lang_B", "Method",
        "Accuracy", "Precision", "Recall", "F1", "Status",
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
