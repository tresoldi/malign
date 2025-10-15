"""Analyze the Arca Verborum gold standard alignment data.

This script analyzes resources/forms_with_alignments.csv to determine:
- Dataset statistics
- Cognate set size distributions
- Suitable subsets for testing
- Data quality metrics
"""

import csv
from collections import Counter, defaultdict
from pathlib import Path


def analyze_gold_data():
    """Analyze the gold standard alignment data."""
    csv_path = Path(__file__).parent.parent / "resources" / "forms_with_alignments.csv"

    datasets = Counter()
    cognate_sets = defaultdict(list)
    cognate_set_sizes = Counter()
    cognate_set_datasets = defaultdict(set)
    total_forms = 0
    alignments_with_gaps = 0

    print(f"Reading {csv_path}...")
    with open(csv_path, encoding="utf-8") as f:
        reader = csv.DictReader(f)

        for row in reader:
            total_forms += 1
            dataset = row["Dataset"]
            cognacy = row["Cognacy"].split(";")[0]  # First cognate set ID
            segments = row["Segments"]
            alignment = row["Alignment"]
            language_id = row["Language_ID"]

            datasets[dataset] += 1
            cognate_set_datasets[cognacy].add(dataset)

            # Store cognate set info
            cognate_sets[cognacy].append({
                "dataset": dataset,
                "language_id": language_id,
                "segments": segments,
                "alignment": alignment,
                "parameter_id": row["Parameter_ID"],
                "gloss": row["Concepticon_Gloss"],
            })

            # Check if alignment has gaps
            if " - " in alignment or alignment.startswith("- ") or alignment.endswith(" -"):
                alignments_with_gaps += 1

    # Calculate cognate set sizes
    for cognacy, forms in cognate_sets.items():
        cognate_set_sizes[len(forms)] += 1

    # Print statistics
    print("\n" + "=" * 80)
    print("GOLD DATA STATISTICS")
    print("=" * 80)

    print(f"\nTotal forms: {total_forms:,}")
    print(f"Total cognate sets: {len(cognate_sets):,}")
    print(f"Total datasets: {len(datasets)}")
    print(f"Forms with gaps in alignment: {alignments_with_gaps:,} ({100*alignments_with_gaps/total_forms:.1f}%)")

    print(f"\n--- Top 10 Datasets by Form Count ---")
    for dataset, count in datasets.most_common(10):
        print(f"  {dataset:30} {count:6,} forms")

    print(f"\n--- Cognate Set Size Distribution ---")
    print(f"  {'Size':>5} {'Count':>10} {'Percentage':>12}")
    for size in sorted(cognate_set_sizes.keys())[:20]:  # Show first 20 sizes
        count = cognate_set_sizes[size]
        pct = 100 * count / len(cognate_sets)
        print(f"  {size:5,} {count:10,} {pct:11.2f}%")

    # Find cognate sets with gaps
    cognate_sets_with_gaps = []
    for cognacy, forms in cognate_sets.items():
        has_gaps = any(" - " in f["alignment"] or f["alignment"].startswith("- ") or f["alignment"].endswith(" -") for f in forms)
        if has_gaps and len(forms) >= 3 and len(forms) <= 10:  # Good size for testing
            cognate_sets_with_gaps.append((cognacy, forms))

    print(f"\n--- Cognate Sets Suitable for Testing (3-10 forms, with gaps) ---")
    print(f"Total suitable sets: {len(cognate_sets_with_gaps):,}")

    # Show examples from different datasets
    print(f"\n--- Example Cognate Sets (first 5 with gaps) ---")
    for i, (cognacy, forms) in enumerate(cognate_sets_with_gaps[:5]):
        if forms:
            print(f"\n  Cognate Set #{i+1}: {cognacy}")
            print(f"    Dataset: {forms[0]['dataset']}")
            print(f"    Gloss: {forms[0]['gloss']}")
            print(f"    Forms: {len(forms)}")
            for j, form in enumerate(forms[:5]):  # Show first 5 forms
                print(f"      {j+1}. {form['language_id']:20} Seg: {form['segments']}")
                print(f"         {' '*20} Aln: {form['alignment']}")
            if len(forms) > 5:
                print(f"      ... and {len(forms)-5} more forms")

    # Recommend datasets for testing
    print(f"\n--- Recommended Datasets for Testing ---")

    # Find datasets with good coverage (many cognate sets of reasonable size)
    dataset_quality = {}
    for dataset in datasets:
        # Count cognate sets in this dataset with 3-10 forms
        good_sets = sum(1 for cognacy, forms in cognate_sets.items()
                       if any(f["dataset"] == dataset for f in forms) and 3 <= len(forms) <= 10)
        dataset_quality[dataset] = good_sets

    for dataset, good_sets in sorted(dataset_quality.items(), key=lambda x: x[1], reverse=True)[:10]:
        print(f"  {dataset:30} {good_sets:6,} suitable cognate sets")

    return cognate_sets, datasets


if __name__ == "__main__":
    analyze_gold_data()
