"""Regression tests using gold standard alignments from Arca Verborum.

These tests verify that MAlign maintains quality on real linguistic data
by comparing predicted alignments against expert-curated gold standards.
"""

from pathlib import Path

import pytest

import malign
from malign.metrics import alignment_accuracy, alignment_f1

from tests.gold_data_utils import (
    cognate_set_to_gold_alignment,
    cognate_set_to_sequences,
    load_cognate_sets,
)


# Path to gold data
DATA_DIR = Path(__file__).parent / "data" / "cognates"
REGRESSION_SET = DATA_DIR / "regression_test_set.yml"


@pytest.mark.regression
def test_regression_smoke_test():
    """Test: Regression smoke test - quick sanity check on small subset.

    Tests a small subset (≤4 forms) of cognate sets to verify basic
    regression testing infrastructure works. This is a fast smoke test.
    """
    # Load gold cognate sets
    cognate_sets = load_cognate_sets(REGRESSION_SET)
    assert len(cognate_sets) > 0, "No cognate sets loaded"

    # Track accuracies
    accuracies = []
    successful_tests = 0
    max_tests = 20  # Test only first 20 that meet criteria

    for cog_set in cognate_sets:
        if successful_tests >= max_tests:
            break

        # Skip empty or malformed cognate sets
        if len(cog_set.forms) < 2 or len(cog_set.forms) > 4:
            continue

        # Skip sets with empty segments/alignments
        if any(not form.segments or not form.alignment for form in cog_set.forms):
            continue

        try:
            # Get sequences and gold alignment
            sequences = cognate_set_to_sequences(cog_set)
            gold_alignment = cognate_set_to_gold_alignment(cog_set)

            # Align using ANW method
            predicted_alms = malign.multi_align(sequences, k=1, method="anw")

            if predicted_alms:
                predicted_alignment = predicted_alms[0]
                acc = alignment_accuracy(predicted_alignment, gold_alignment)
                accuracies.append(acc)
                successful_tests += 1

        except Exception:
            continue

    # Verify we tested a reasonable number of sets
    assert successful_tests >= 15, f"Only {successful_tests} tests succeeded (expected >= 15)"

    # Compute average metrics
    avg_accuracy = sum(accuracies) / len(accuracies)

    # Report results
    print(f"\nRegression smoke test: {successful_tests} sets, avg accuracy: {avg_accuracy:.2%}")

    # Assert 60% threshold (baseline - room for improvement with learned matrices)
    # Note: Default identity matrix achieves ~63% on this gold data.
    # Future improvements: matrix learning from gold alignments should improve this.
    assert avg_accuracy >= 0.60, (
        f"Average accuracy {avg_accuracy:.2%} below 60% baseline threshold. "
        f"Tested {successful_tests} cognate sets."
    )


@pytest.mark.regression
@pytest.mark.slow
def test_regression_accuracy_threshold_anw():
    """Test: Regression - ANW method baseline accuracy on gold data.

    This test loads cognate sets from the Arca Verborum gold standard,
    aligns them using the ANW (A-star + Needleman-Wunsch) method with
    default identity matrix, and verifies baseline accuracy (65%+).

    Note: 65% baseline with identity matrix. Matrix learning from gold data
    should improve this significantly. Marked @pytest.mark.slow as computationally intensive.
    """
    # Load gold cognate sets
    cognate_sets = load_cognate_sets(REGRESSION_SET)
    assert len(cognate_sets) > 0, "No cognate sets loaded"

    # Track accuracies
    accuracies = []
    f1_scores = []
    successful_tests = 0
    failed_tests = []

    for cog_set in cognate_sets:
        # Skip empty or malformed cognate sets
        if len(cog_set.forms) < 2:
            continue

        # Skip large multiwise sets (> 5 forms) for reasonable runtime
        if len(cog_set.forms) > 5:
            continue

        # Skip sets with empty segments/alignments
        if any(not form.segments or not form.alignment for form in cog_set.forms):
            continue

        try:
            # Get sequences and gold alignment
            sequences = cognate_set_to_sequences(cog_set)
            gold_alignment = cognate_set_to_gold_alignment(cog_set)

            # Align using ANW method
            predicted_alms = malign.multi_align(sequences, k=1, method="anw")

            if not predicted_alms:
                failed_tests.append((cog_set.id, "No alignments produced"))
                continue

            predicted_alignment = predicted_alms[0]

            # Compute metrics
            acc = alignment_accuracy(predicted_alignment, gold_alignment)
            f1 = alignment_f1(predicted_alignment, gold_alignment)

            accuracies.append(acc)
            f1_scores.append(f1)
            successful_tests += 1

        except Exception as e:
            # Record failures but continue
            failed_tests.append((cog_set.id, str(e)))

    # Verify we tested a reasonable number of sets
    assert successful_tests >= 60, f"Only {successful_tests} tests succeeded (expected >= 60)"

    # Compute average metrics
    avg_accuracy = sum(accuracies) / len(accuracies)
    avg_f1 = sum(f1_scores) / len(f1_scores)

    # Report results
    print(f"\n{'='*80}")
    print(f"REGRESSION TEST RESULTS (ANW Method)")
    print(f"{'='*80}")
    print(f"Total cognate sets: {len(cognate_sets)}")
    print(f"Successfully tested: {successful_tests}")
    print(f"Failed tests: {len(failed_tests)}")
    print(f"Average accuracy: {avg_accuracy:.2%}")
    print(f"Average F1 score: {avg_f1:.2%}")
    print(f"Accuracy range: [{min(accuracies):.2%}, {max(accuracies):.2%}]")
    print(f"{'='*80}")

    if failed_tests:
        print(f"\nFailed tests (first 10):")
        for cog_id, error in failed_tests[:10]:
            print(f"  {cog_id}: {error}")

    # Assert 65% baseline threshold
    assert avg_accuracy >= 0.65, (
        f"Average accuracy {avg_accuracy:.2%} below 65% baseline threshold. "
        f"Tested {successful_tests} cognate sets. "
        f"Note: Identity matrix baseline. Learned matrices should improve this."
    )


@pytest.mark.regression
@pytest.mark.slow
def test_regression_accuracy_threshold_yenksp():
    """Test: Regression - YenKSP method baseline accuracy on gold data.

    This test verifies baseline accuracy (65%+) using the YenKSP
    (Yen's k-shortest paths) method with default identity matrix.

    Note: 65% baseline. Matrix learning should improve this.
    Marked @pytest.mark.slow as YenKSP is more computationally intensive.
    """
    # Load gold cognate sets
    cognate_sets = load_cognate_sets(REGRESSION_SET)
    assert len(cognate_sets) > 0, "No cognate sets loaded"

    # Track accuracies
    accuracies = []
    f1_scores = []
    successful_tests = 0
    failed_tests = []

    for cog_set in cognate_sets:
        # Skip empty or malformed cognate sets
        if len(cog_set.forms) < 2:
            continue

        # Skip sets with empty segments/alignments
        if any(not form.segments or not form.alignment for form in cog_set.forms):
            continue

        # Skip sets with > 4 forms (YenKSP is slow on multiwise)
        if len(cog_set.forms) > 4:
            continue

        try:
            # Get sequences and gold alignment
            sequences = cognate_set_to_sequences(cog_set)
            gold_alignment = cognate_set_to_gold_alignment(cog_set)

            # Align using YenKSP method
            predicted_alms = malign.multi_align(sequences, k=1, method="yenksp")

            if not predicted_alms:
                failed_tests.append((cog_set.id, "No alignments produced"))
                continue

            predicted_alignment = predicted_alms[0]

            # Compute metrics
            acc = alignment_accuracy(predicted_alignment, gold_alignment)
            f1 = alignment_f1(predicted_alignment, gold_alignment)

            accuracies.append(acc)
            f1_scores.append(f1)
            successful_tests += 1

        except Exception as e:
            # Record failures but continue
            failed_tests.append((cog_set.id, str(e)))

    # Verify we tested a reasonable number of sets
    assert successful_tests >= 40, f"Only {successful_tests} tests succeeded (expected >= 40)"

    # Compute average metrics
    avg_accuracy = sum(accuracies) / len(accuracies)
    avg_f1 = sum(f1_scores) / len(f1_scores)

    # Report results
    print(f"\n{'='*80}")
    print(f"REGRESSION TEST RESULTS (YenKSP Method)")
    print(f"{'='*80}")
    print(f"Total cognate sets: {len(cognate_sets)}")
    print(f"Successfully tested: {successful_tests}")
    print(f"Failed tests: {len(failed_tests)}")
    print(f"Average accuracy: {avg_accuracy:.2%}")
    print(f"Average F1 score: {avg_f1:.2%}")
    print(f"Accuracy range: [{min(accuracies):.2%}, {max(accuracies):.2%}]")
    print(f"{'='*80}")

    if failed_tests:
        print(f"\nFailed tests (first 10):")
        for cog_id, error in failed_tests[:10]:
            print(f"  {cog_id}: {error}")

    # Assert 65% baseline threshold
    assert avg_accuracy >= 0.65, (
        f"Average accuracy {avg_accuracy:.2%} below 65% baseline threshold. "
        f"Tested {successful_tests} cognate sets."
    )


@pytest.mark.regression
def test_regression_per_dataset_quality():
    """Test: Regression - Verify quality across different datasets.

    This test analyzes alignment quality broken down by dataset,
    ensuring consistent performance across different linguistic families.
    """
    # Load gold cognate sets
    cognate_sets = load_cognate_sets(REGRESSION_SET)

    # Group by dataset
    by_dataset = {}
    for cog_set in cognate_sets:
        if cog_set.dataset not in by_dataset:
            by_dataset[cog_set.dataset] = []
        by_dataset[cog_set.dataset].append(cog_set)

    dataset_results = {}

    for dataset, sets in by_dataset.items():
        accuracies = []

        for cog_set in sets:
            if len(cog_set.forms) < 2 or len(cog_set.forms) > 5:
                continue

            if any(not form.segments or not form.alignment for form in cog_set.forms):
                continue

            try:
                sequences = cognate_set_to_sequences(cog_set)
                gold_alignment = cognate_set_to_gold_alignment(cog_set)
                predicted_alms = malign.multi_align(sequences, k=1, method="anw")

                if predicted_alms:
                    acc = alignment_accuracy(predicted_alms[0], gold_alignment)
                    accuracies.append(acc)
            except Exception:
                continue

        if accuracies:
            avg_acc = sum(accuracies) / len(accuracies)
            dataset_results[dataset] = {
                "count": len(accuracies),
                "avg_accuracy": avg_acc,
            }

    # Report per-dataset results
    print(f"\n{'='*80}")
    print(f"PER-DATASET ACCURACY")
    print(f"{'='*80}")
    for dataset, results in sorted(dataset_results.items(), key=lambda x: x[1]["avg_accuracy"], reverse=True):
        print(f"{dataset:30} {results['count']:3} sets, {results['avg_accuracy']:.2%} accuracy")
    print(f"{'='*80}")

    # Verify at least 3 datasets tested
    assert len(dataset_results) >= 3, f"Only {len(dataset_results)} datasets tested"

    # Verify each dataset has reasonable quality (>= 50% baseline)
    # Note: Identity matrix baseline. Variation across linguistic families is expected.
    for dataset, results in dataset_results.items():
        assert results["avg_accuracy"] >= 0.50, (
            f"Dataset {dataset} has very low accuracy: {results['avg_accuracy']:.2%}"
        )
