"""Tests for alignment validation metrics (Phase 2.3)."""

import pytest

import malign
from malign.alignment import Alignment


def test_alignment_accuracy_perfect_match():
    """Test accuracy with perfectly matching alignments."""
    predicted = Alignment([["A", "C", "G", "T"], ["A", "C", "G", "T"]], score=1.0)
    gold = Alignment([["A", "C", "G", "T"], ["A", "C", "G", "T"]], score=1.0)

    accuracy = malign.alignment_accuracy(predicted, gold)
    assert accuracy == 1.0


def test_alignment_accuracy_partial_match():
    """Test accuracy with partially matching alignments."""
    # Columns: [A,A], [C,C], [-,G], [T,T]
    # Matches: column 0, 1, 3 = 3/4 = 0.75
    predicted = Alignment([["A", "C", "-", "T"], ["A", "C", "G", "T"]], score=1.0)
    gold = Alignment([["A", "C", "G", "T"], ["A", "C", "G", "T"]], score=1.0)

    accuracy = malign.alignment_accuracy(predicted, gold)
    assert accuracy == 0.75


def test_alignment_accuracy_no_match():
    """Test accuracy with completely mismatched alignments."""
    predicted = Alignment([["A", "A"], ["C", "C"]], score=1.0)
    gold = Alignment([["C", "C"], ["A", "A"]], score=1.0)

    accuracy = malign.alignment_accuracy(predicted, gold)
    assert accuracy == 0.0


def test_alignment_accuracy_empty():
    """Test accuracy with empty alignments."""
    predicted = Alignment([[], []], score=0.0)
    gold = Alignment([[], []], score=0.0)

    accuracy = malign.alignment_accuracy(predicted, gold)
    assert accuracy == 1.0


def test_alignment_accuracy_length_mismatch():
    """Test that length mismatch raises ValueError."""
    predicted = Alignment([["A", "C"], ["A", "C"]], score=1.0)
    gold = Alignment([["A", "C", "G"], ["A", "C", "G"]], score=1.0)

    with pytest.raises(ValueError, match="length mismatch"):
        malign.alignment_accuracy(predicted, gold)


def test_alignment_accuracy_sequence_count_mismatch():
    """Test that sequence count mismatch raises ValueError."""
    predicted = Alignment([["A", "C"], ["A", "C"]], score=1.0)
    gold = Alignment([["A", "C"], ["A", "C"], ["A", "C"]], score=1.0)

    with pytest.raises(ValueError, match="Sequence count mismatch"):
        malign.alignment_accuracy(predicted, gold)


def test_precision_recall_perfect_match():
    """Test precision/recall with perfect match."""
    predicted = Alignment([["A", "C"], ["A", "C"]], score=1.0)
    gold = Alignment([["A", "C"], ["A", "C"]], score=1.0)

    precision, recall = malign.alignment_precision_recall(predicted, gold)
    assert precision == 1.0
    assert recall == 1.0


def test_precision_recall_partial_match():
    """Test precision/recall with partial match."""
    # Predicted has gap in first sequence, gold doesn't
    predicted = Alignment([["-", "C"], ["A", "C"]], score=1.0)
    gold = Alignment([["A", "C"], ["A", "C"]], score=1.0)

    precision, recall = malign.alignment_precision_recall(predicted, gold)

    # Not all predicted pairs are in gold (lower precision)
    # Not all gold pairs are predicted (lower recall)
    assert 0.0 < precision < 1.0
    assert 0.0 < recall < 1.0


def test_precision_recall_no_overlap():
    """Test precision/recall with no overlap."""
    predicted = Alignment([["A", "A"], ["C", "C"]], score=1.0)
    gold = Alignment([["C", "C"], ["A", "A"]], score=1.0)

    precision, recall = malign.alignment_precision_recall(predicted, gold)

    # Completely different alignments = no overlap
    assert precision == 0.0
    assert recall == 0.0


def test_precision_recall_empty():
    """Test precision/recall with empty alignments."""
    predicted = Alignment([[], []], score=0.0)
    gold = Alignment([[], []], score=0.0)

    precision, recall = malign.alignment_precision_recall(predicted, gold)
    assert precision == 1.0
    assert recall == 1.0


def test_f1_perfect_match():
    """Test F1 score with perfect match."""
    predicted = Alignment([["A", "C", "G"], ["A", "C", "G"]], score=1.0)
    gold = Alignment([["A", "C", "G"], ["A", "C", "G"]], score=1.0)

    f1 = malign.alignment_f1(predicted, gold)
    assert f1 == 1.0


def test_f1_partial_match():
    """Test F1 score with partial match."""
    predicted = Alignment([["-", "C"], ["A", "C"]], score=1.0)
    gold = Alignment([["A", "C"], ["A", "C"]], score=1.0)

    f1 = malign.alignment_f1(predicted, gold)

    # F1 should be between 0 and 1, and is harmonic mean of P and R
    assert 0.0 < f1 < 1.0


def test_f1_no_match():
    """Test F1 score with no match."""
    predicted = Alignment([["A", "A"], ["C", "C"]], score=1.0)
    gold = Alignment([["C", "C"], ["A", "A"]], score=1.0)

    f1 = malign.alignment_f1(predicted, gold)
    assert f1 == 0.0


def test_f1_is_harmonic_mean():
    """Test that F1 is correctly computed as harmonic mean."""
    predicted = Alignment([["A", "C", "G"], ["A", "-", "G"]], score=1.0)
    gold = Alignment([["A", "C", "G"], ["A", "C", "G"]], score=1.0)

    precision, recall = malign.alignment_precision_recall(predicted, gold)
    f1 = malign.alignment_f1(predicted, gold)

    # Verify F1 = 2 * (P * R) / (P + R)
    expected_f1 = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0.0
    assert abs(f1 - expected_f1) < 1e-10


def test_metrics_with_multisequence_alignment():
    """Test metrics with more than 2 sequences."""
    predicted = Alignment([["A", "C"], ["A", "C"], ["A", "C"]], score=1.0)
    gold = Alignment([["A", "C"], ["A", "C"], ["A", "C"]], score=1.0)

    accuracy = malign.alignment_accuracy(predicted, gold)
    precision, recall = malign.alignment_precision_recall(predicted, gold)
    f1 = malign.alignment_f1(predicted, gold)

    assert accuracy == 1.0
    assert precision == 1.0
    assert recall == 1.0
    assert f1 == 1.0
