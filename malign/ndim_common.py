"""Common infrastructure for N-dimensional alignment."""

import math
from collections.abc import Hashable

# Thresholds for true N-dim alignment feasibility.
# N-dim always uses YenKSP (NW backtrace has exponential co-optimal explosion).
_NDIM_MAX_N = 4
_NDIM_MAX_CELLS = 200_000


def make_moves(n: int) -> list[tuple[int, tuple[int, ...]]]:
    """Generate all 2^N - 1 possible moves in an N-dimensional grid.

    Each move is a (bitmask, delta_tuple) pair. Bit d is set if dimension d
    advances; delta[d] is -1 for advancing dims, 0 otherwise.

    Args:
        n: Number of dimensions.

    Returns:
        List of (bitmask, delta) pairs for all non-zero moves.
    """
    moves = []
    for bitmask in range(1, 1 << n):
        delta = tuple(-1 if bitmask & (1 << d) else 0 for d in range(n))
        moves.append((bitmask, delta))
    return moves


def should_use_ndim(
    n: int,
    seq_lengths: list[int],
    method: str,  # noqa: ARG001
) -> bool:
    """Check whether true N-dimensional alignment is feasible.

    N-dim alignment always uses YenKSP internally (Dijkstra-based,
    produces exactly k paths without backtrace explosion).

    Args:
        n: Number of sequences.
        seq_lengths: Length of each sequence.
        method: Alignment method (unused, kept for API compatibility).

    Returns:
        True if N-dim alignment is feasible, False for progressive fallback.
    """
    grid_cells = math.prod(length + 1 for length in seq_lengths)
    return n <= _NDIM_MAX_N and grid_cells <= _NDIM_MAX_CELLS


def _column_symbols(
    seqs: list[list[Hashable]],
    idx: tuple[int, ...],
    bitmask: int,
    gap: Hashable,
) -> tuple[Hashable, ...]:
    """Build the symbol tuple for a move at a given grid position.

    Args:
        seqs: Padded sequences (with leading gap).
        idx: Current N-dimensional grid index.
        bitmask: Move bitmask indicating which dims advance.
        gap: Gap symbol.

    Returns:
        Tuple of symbols for the alignment column.
    """
    n = len(seqs)
    return tuple(seqs[d][idx[d]] if bitmask & (1 << d) else gap for d in range(n))
