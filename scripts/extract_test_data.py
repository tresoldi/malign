"""Extract test data subsets from Arca Verborum gold standard.

This script extracts curated subsets for:
- Regression tests (50-100 cognate sets with gold alignments)
- Matrix learning training (100-200 cognate sets)
- Matrix learning evaluation (50-100 different sets)
- Integration tests (10-20 diverse examples)
"""

import csv
import random
import yaml
from collections import defaultdict
from pathlib import Path


def load_cognate_sets(csv_path):
    """Load all cognate sets from CSV."""
    cognate_sets = defaultdict(list)

    with open(csv_path, encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            cognacy = row["Cognacy"].split(";")[0]
            cognate_sets[cognacy].append({
                "dataset": row["Dataset"],
                "language_id": row["Language_ID"],
                "glottolog_name": row["Glottolog_Name"],
                "parameter_id": row["Parameter_ID"],
                "gloss": row["Concepticon_Gloss"],
                "segments": row["Segments"].split(),
                "alignment": row["Alignment"].split(),
            })

    return cognate_sets


def has_gaps(cognate_set):
    """Check if cognate set has gaps in alignments."""
    return any("-" in form["alignment"] for form in cognate_set)


def extract_regression_set(cognate_sets, n=100, min_size=3, max_size=10):
    """Extract diverse cognate sets for regression testing."""
    # Filter: size 3-10, has gaps, high quality
    candidates = [
        (cognacy, forms) for cognacy, forms in cognate_sets.items()
        if min_size <= len(forms) <= max_size and has_gaps(forms)
    ]

    # Sample diverse datasets
    by_dataset = defaultdict(list)
    for cognacy, forms in candidates:
        by_dataset[forms[0]["dataset"]].append((cognacy, forms))

    # Take proportionally from each dataset
    selected = []
    datasets = list(by_dataset.keys())
    random.seed(42)  # Reproducible
    random.shuffle(datasets)

    for dataset in datasets:
        items = by_dataset[dataset]
        random.shuffle(items)
        # Take up to 10 from each dataset
        selected.extend(items[:min(10, len(items))])
        if len(selected) >= n:
            break

    return selected[:n]


def extract_learning_training_set(cognate_sets, n=200, min_size=2, max_size=8):
    """Extract cognate sets for matrix learning training."""
    candidates = [
        (cognacy, forms) for cognacy, forms in cognate_sets.items()
        if min_size <= len(forms) <= max_size
    ]

    random.seed(43)
    random.shuffle(candidates)
    return candidates[:n]


def extract_learning_eval_set(cognate_sets, training_set_ids, n=100):
    """Extract different cognate sets for evaluation (not in training)."""
    training_ids = set(cog_id for cog_id, _ in training_set_ids)

    candidates = [
        (cognacy, forms) for cognacy, forms in cognate_sets.items()
        if cognacy not in training_ids and 2 <= len(forms) <= 8 and has_gaps(forms)
    ]

    random.seed(44)
    random.shuffle(candidates)
    return candidates[:n]


def extract_integration_examples(cognate_sets, n=20):
    """Extract diverse examples for integration tests."""
    # Select examples with various characteristics
    examples = []

    # 1. Small sets (2-3 forms)
    small = [(c, f) for c, f in cognate_sets.items() if 2 <= len(f) <= 3 and has_gaps(f)]
    random.seed(45)
    random.shuffle(small)
    examples.extend(small[:5])

    # 2. Medium sets (4-6 forms)
    medium = [(c, f) for c, f in cognate_sets.items() if 4 <= len(f) <= 6 and has_gaps(f)]
    random.shuffle(medium)
    examples.extend(medium[:8])

    # 3. Larger sets (7-10 forms)
    large = [(c, f) for c, f in cognate_sets.items() if 7 <= len(f) <= 10 and has_gaps(f)]
    random.shuffle(large)
    examples.extend(large[:7])

    return examples[:n]


def save_as_yaml(cognate_sets_list, output_path, description):
    """Save cognate sets as YAML."""
    data = {
        "description": description,
        "source": "Arca Verborum (forms_with_alignments.csv)",
        "count": len(cognate_sets_list),
        "cognate_sets": []
    }

    for cognacy_id, forms in cognate_sets_list:
        cognate_set = {
            "id": cognacy_id,
            "dataset": forms[0]["dataset"],
            "gloss": forms[0]["gloss"],
            "parameter_id": forms[0]["parameter_id"],
            "forms": []
        }

        for form in forms:
            cognate_set["forms"].append({
                "language_id": form["language_id"],
                "glottolog_name": form["glottolog_name"],
                "segments": form["segments"],
                "alignment": form["alignment"],
            })

        data["cognate_sets"].append(cognate_set)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as f:
        yaml.dump(data, f, default_flow_style=False, allow_unicode=True, sort_keys=False)

    print(f"Saved {len(cognate_sets_list)} cognate sets to {output_path}")


def main():
    """Extract all test data subsets."""
    csv_path = Path(__file__).parent.parent / "resources" / "forms_with_alignments.csv"
    output_dir = Path(__file__).parent.parent / "tests" / "data" / "cognates"

    print("Loading cognate sets from CSV...")
    cognate_sets = load_cognate_sets(csv_path)
    print(f"Loaded {len(cognate_sets):,} cognate sets")

    print("\nExtracting regression test set (100 cognate sets)...")
    regression_set = extract_regression_set(cognate_sets, n=100)
    save_as_yaml(
        regression_set,
        output_dir / "regression_test_set.yml",
        "Gold standard alignments for regression testing (100 diverse cognate sets)"
    )

    print("\nExtracting matrix learning training set (200 cognate sets)...")
    training_set = extract_learning_training_set(cognate_sets, n=200)
    save_as_yaml(
        training_set,
        output_dir / "learning_training_set.yml",
        "Cognate sets for training matrix learning algorithms (200 sets)"
    )

    print("\nExtracting matrix learning evaluation set (100 cognate sets)...")
    eval_set = extract_learning_eval_set(cognate_sets, training_set, n=100)
    save_as_yaml(
        eval_set,
        output_dir / "learning_eval_set.yml",
        "Cognate sets for evaluating learned matrices (100 different sets)"
    )

    print("\nExtracting integration test examples (20 cognate sets)...")
    integration_examples = extract_integration_examples(cognate_sets, n=20)
    save_as_yaml(
        integration_examples,
        output_dir / "integration_examples.yml",
        "Diverse examples for integration tests (20 sets of varying sizes)"
    )

    print("\n" + "="*80)
    print("EXTRACTION COMPLETE")
    print("="*80)
    print(f"\nFiles created in {output_dir}:")
    print(f"  - regression_test_set.yml (100 cognate sets)")
    print(f"  - learning_training_set.yml (200 cognate sets)")
    print(f"  - learning_eval_set.yml (100 cognate sets)")
    print(f"  - integration_examples.yml (20 cognate sets)")


if __name__ == "__main__":
    main()
