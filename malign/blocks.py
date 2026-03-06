"""Block detection and column merging for phonological alignments.

Detects complementary-gap patterns (diphthongization, metathesis) in
alignments and merges block columns into compound symbols (tuples).
"""

from collections import Counter
from collections.abc import Hashable, Sequence
from enum import Enum, auto

from .alignment import Alignment


class _ColType(Enum):
    FULL = auto()
    PARTIAL = auto()
    ALL_GAP = auto()


def _classify_column(column: Sequence[Hashable], gap: Hashable) -> _ColType:
    """Classify a single alignment column."""
    has_gap = any(s == gap for s in column)
    all_gap = all(s == gap for s in column)
    if all_gap:
        return _ColType.ALL_GAP
    if has_gap:
        return _ColType.PARTIAL
    return _ColType.FULL


def detect_blocks(
    seqs: Sequence[Sequence[Hashable]],
    gap: Hashable = "-",
    max_block_size: int = 2,
) -> list[tuple[int, int]]:
    """Detect block patterns (complementary-gap columns) in an alignment.

    Finds maximal runs of consecutive PARTIAL columns and extends them
    to adjacent FULL columns (one per side) when possible.

    Args:
        seqs: Aligned sequences (all same length).
        gap: Gap symbol (default: "-").
        max_block_size: Maximum number of columns in a block (default: 2).

    Returns:
        List of (start, end) inclusive index tuples for each detected block.
    """
    if not seqs or not seqs[0]:
        return []

    num_cols = len(seqs[0])
    col_types = [_classify_column(tuple(seq[i] for seq in seqs), gap) for i in range(num_cols)]

    blocks: list[tuple[int, int]] = []
    i = 0
    while i < num_cols:
        # Skip non-PARTIAL columns
        if col_types[i] != _ColType.PARTIAL:
            i += 1
            continue

        # Find maximal run of consecutive PARTIALs
        run_start = i
        while i < num_cols and col_types[i] == _ColType.PARTIAL:
            i += 1
        run_end = i - 1  # inclusive

        run_len = run_end - run_start + 1

        # If the run itself already fits within max_block_size, try extending
        if run_len <= max_block_size:
            can_extend_left = run_start > 0 and col_types[run_start - 1] == _ColType.FULL
            can_extend_right = run_end < num_cols - 1 and col_types[run_end + 1] == _ColType.FULL

            if can_extend_left and can_extend_right:
                # Both sides possible — try both first
                if run_len + 2 <= max_block_size:
                    blocks.append((run_start - 1, run_end + 1))
                else:
                    # Can't extend both — extend left only
                    if run_len + 1 <= max_block_size:
                        blocks.append((run_start - 1, run_end))
                    else:
                        # Run itself is the block
                        if run_len >= 2:
                            blocks.append((run_start, run_end))
            elif can_extend_left:
                if run_len + 1 <= max_block_size:
                    blocks.append((run_start - 1, run_end))
                elif run_len >= 2:
                    blocks.append((run_start, run_end))
            elif can_extend_right:
                if run_len + 1 <= max_block_size:
                    blocks.append((run_start, run_end + 1))
                elif run_len >= 2:
                    blocks.append((run_start, run_end))
            else:
                # No adjacent FULL — block only if run itself is >= 2
                if run_len >= 2:
                    blocks.append((run_start, run_end))
        elif run_len > max_block_size:
            # Run too long — no valid block
            pass

    return blocks


def merge_alignment_blocks(
    alignment: Alignment,
    gap: Hashable = "-",
    max_block_size: int = 2,
) -> Alignment:
    """Merge block columns into compound symbols.

    Non-gap symbols within a block become tuples; single symbols stay
    unwrapped; all-gap positions become the gap symbol.

    Args:
        alignment: The alignment to process.
        gap: Gap symbol (default: "-").
        max_block_size: Maximum block size for detection (default: 2).

    Returns:
        New Alignment with merged block columns.
    """
    seqs = alignment.seqs
    if not seqs or not seqs[0]:
        return alignment

    blocks = detect_blocks(seqs, gap=gap, max_block_size=max_block_size)

    # Build set of column indices that belong to blocks
    block_cols: set[int] = set()
    for start, end in blocks:
        block_cols.update(range(start, end + 1))

    num_seqs = len(seqs)
    num_cols = len(seqs[0])

    # Build new columns
    new_columns: list[tuple[Hashable, ...]] = []
    i = 0
    while i < num_cols:
        if i not in block_cols:
            # Non-block column — pass through
            new_columns.append(tuple(seq[i] for seq in seqs))
            i += 1
        else:
            # Find which block this belongs to
            for start, end in blocks:
                if start <= i <= end:
                    # Merge all columns in this block
                    col = []
                    for seq_idx in range(num_seqs):
                        symbols = [
                            seqs[seq_idx][c]
                            for c in range(start, end + 1)
                            if seqs[seq_idx][c] != gap
                        ]
                        if len(symbols) == 0:
                            col.append(gap)
                        elif len(symbols) == 1:
                            col.append(symbols[0])
                        else:
                            col.append(tuple(symbols))
                    new_columns.append(tuple(col))
                    i = end + 1
                    break

    # Reconstruct sequences from columns
    new_seqs = [tuple(col[seq_idx] for col in new_columns) for seq_idx in range(num_seqs)]

    return Alignment(seqs=new_seqs, score=alignment.score)


def _adjust_pair_counts_for_blocks(
    pair_counts: Counter[tuple[Hashable, ...]],
    alignment: Alignment,
    gap: Hashable,
    max_block_size: int,
) -> None:
    """Decrement pair_counts for gap-containing columns inside blocks.

    Modifies pair_counts in place. Only PARTIAL columns within blocks
    are decremented; FULL columns are left untouched. Counts are clamped
    to zero.

    Args:
        pair_counts: Counter of column tuples to adjust.
        alignment: The alignment whose blocks to analyze.
        gap: Gap symbol.
        max_block_size: Maximum block size for detection.
    """
    seqs = alignment.seqs
    if not seqs or not seqs[0]:
        return

    blocks = detect_blocks(seqs, gap=gap, max_block_size=max_block_size)

    for start, end in blocks:
        for col_idx in range(start, end + 1):
            column = tuple(seq[col_idx] for seq in seqs)
            col_type = _classify_column(column, gap)
            if col_type == _ColType.PARTIAL:
                if pair_counts[column] > 0:
                    pair_counts[column] -= 1
                # Clean up zero entries
                if pair_counts[column] <= 0:
                    pair_counts.pop(column, None)
