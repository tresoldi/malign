"""Validation metrics for alignment quality assessment.

This module provides metrics to evaluate alignment quality by comparing
predicted alignments against gold-standard reference alignments.
"""

from collections.abc import Hashable

from .alignment import Alignment


def strip_common_gaps(alignment: Alignment, gap: Hashable = "-") -> Alignment:
    """Remove columns where all sequences are gaps.

    When gold alignments are extracted from a multiple sequence alignment (MSA),
    pairwise slices may contain columns where both sequences are gaps (induced
    by a third sequence in the MSA). No pairwise aligner will produce these
    all-gap columns, so they must be stripped before comparing a predicted
    pairwise alignment against an MSA-derived gold standard.

    Args:
        alignment: The alignment to strip.
        gap: Gap symbol (default: "-").

    Returns:
        A new Alignment with all-gap columns removed and the same score.

    Example:
        >>> gold = Alignment([["a", "-", "-", "b"], ["a", "-", "c", "b"]], score=0.0)
        >>> stripped = strip_common_gaps(gold)
        >>> stripped.seqs  # Column 1 (all gaps) removed; column 2 kept (only one gap)
        [['a', '-', 'b'], ['a', 'c', 'b']]
    """
    if not alignment.seqs or not alignment.seqs[0]:
        return alignment

    num_seqs = len(alignment.seqs)
    aln_len = len(alignment.seqs[0])

    new_seqs: list[list[Hashable]] = [[] for _ in range(num_seqs)]
    for col_idx in range(aln_len):
        column = [alignment.seqs[seq_idx][col_idx] for seq_idx in range(num_seqs)]
        if all(s == gap for s in column):
            continue
        for seq_idx in range(num_seqs):
            new_seqs[seq_idx].append(column[seq_idx])

    return Alignment(seqs=new_seqs, score=alignment.score)


def alignment_accuracy(predicted: Alignment, gold: Alignment) -> float:
    """Calculate alignment accuracy as the proportion of matching columns.

    Compares two alignments column-by-column. A column matches if all sequences
    at that position have identical symbols (including gaps).

    Args:
        predicted: The predicted alignment to evaluate.
        gold: The gold-standard reference alignment.

    Returns:
        Accuracy score in [0.0, 1.0], where 1.0 means perfect match.

    Raises:
        ValueError: If alignments have different lengths or sequence counts.

    Example:
        >>> predicted = Alignment([["A", "C", "-"], ["A", "C", "T"]], score=1.0)
        >>> gold = Alignment([["A", "C", "G"], ["A", "C", "T"]], score=1.0)
        >>> accuracy = alignment_accuracy(predicted, gold)
        >>> # Columns 0 and 1 match, column 2 doesn't: accuracy = 2/3 ≈ 0.667
    """
    # Validate inputs
    if len(predicted.seqs) != len(gold.seqs):
        raise ValueError(
            f"Sequence count mismatch: predicted has {len(predicted.seqs)}, "
            f"gold has {len(gold.seqs)}"
        )

    # Get alignment lengths
    pred_len = len(predicted.seqs[0]) if predicted.seqs else 0
    gold_len = len(gold.seqs[0]) if gold.seqs else 0

    if pred_len != gold_len:
        raise ValueError(
            f"Alignment length mismatch: predicted has {pred_len}, gold has {gold_len}"
        )

    if pred_len == 0:
        return 1.0  # Empty alignments match perfectly

    # Compare column by column
    matching_columns = 0
    for col_idx in range(pred_len):
        pred_col = tuple(seq[col_idx] for seq in predicted.seqs)
        gold_col = tuple(seq[col_idx] for seq in gold.seqs)
        if pred_col == gold_col:
            matching_columns += 1

    return matching_columns / pred_len


def _get_alignment_pairs(alignment: Alignment) -> set[tuple[int, int, Hashable, Hashable]]:
    """Extract aligned symbol pairs from an alignment.

    For each column, extracts all pairwise symbol alignments between sequences.
    This is used for precision/recall calculations.

    Args:
        alignment: The alignment to extract pairs from.

    Returns:
        Set of (seq_i, seq_j, symbol_i, symbol_j, col_idx) tuples representing
        all pairwise alignments in the alignment.
    """
    pairs = set()
    num_seqs = len(alignment.seqs)
    aln_len = len(alignment.seqs[0]) if alignment.seqs else 0

    for col_idx in range(aln_len):
        # Get all symbols at this column
        column = [seq[col_idx] for seq in alignment.seqs]

        # Extract all pairwise alignments from this column
        for i in range(num_seqs):
            for j in range(i + 1, num_seqs):
                pairs.add((i, j, column[i], column[j], col_idx))

    return pairs


def alignment_precision_recall(predicted: Alignment, gold: Alignment) -> tuple[float, float]:
    """Calculate precision and recall for an alignment.

    Precision: proportion of predicted pairwise alignments that are correct.
    Recall: proportion of gold pairwise alignments that were predicted.

    This treats alignments as sets of pairwise symbol alignments and compares
    them using standard information retrieval metrics.

    Args:
        predicted: The predicted alignment to evaluate.
        gold: The gold-standard reference alignment.

    Returns:
        Tuple of (precision, recall), both in [0.0, 1.0].

    Raises:
        ValueError: If alignments have different sequence counts.

    Example:
        >>> predicted = Alignment([["A", "C"], ["A", "C"]], score=1.0)
        >>> gold = Alignment([["A", "C"], ["A", "C"]], score=1.0)
        >>> precision, recall = alignment_precision_recall(predicted, gold)
        >>> # Perfect match: precision = recall = 1.0
    """
    # Validate inputs
    if len(predicted.seqs) != len(gold.seqs):
        raise ValueError(
            f"Sequence count mismatch: predicted has {len(predicted.seqs)}, "
            f"gold has {len(gold.seqs)}"
        )

    # Extract pairwise alignments
    pred_pairs = _get_alignment_pairs(predicted)
    gold_pairs = _get_alignment_pairs(gold)

    # Handle edge cases
    if not pred_pairs and not gold_pairs:
        return 1.0, 1.0  # Both empty = perfect match

    if not pred_pairs:
        return 0.0, 0.0  # Predicted nothing

    if not gold_pairs:
        return 0.0, 0.0  # Gold is empty (shouldn't happen in practice)

    # Calculate metrics
    correct = len(pred_pairs & gold_pairs)
    precision = correct / len(pred_pairs)
    recall = correct / len(gold_pairs)

    return precision, recall


def alignment_f1(predicted: Alignment, gold: Alignment) -> float:
    """Calculate F1 score for an alignment.

    F1 is the harmonic mean of precision and recall, providing a single
    metric that balances both measures.

    Args:
        predicted: The predicted alignment to evaluate.
        gold: The gold-standard reference alignment.

    Returns:
        F1 score in [0.0, 1.0], where 1.0 means perfect match.

    Raises:
        ValueError: If alignments have different sequence counts.

    Example:
        >>> predicted = Alignment([["A", "-"], ["A", "C"]], score=1.0)
        >>> gold = Alignment([["A", "C"], ["A", "C"]], score=1.0)
        >>> f1 = alignment_f1(predicted, gold)
    """
    precision, recall = alignment_precision_recall(predicted, gold)

    # Handle edge case where both are zero
    if precision + recall == 0:
        return 0.0

    # Calculate harmonic mean
    f1 = 2 * (precision * recall) / (precision + recall)
    return f1
