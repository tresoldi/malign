"""N-dimensional Needleman-Wunsch alignment."""

from collections.abc import Hashable

import numpy as np

from .alignment import Alignment
from .ndim_common import _column_symbols, make_moves
from .scoring_matrix import ScoringMatrix
from .utils import score_alignment, sort_alignments

_EPS = 1e-12


def ndim_nw_grids(
    seqs: list[list[Hashable]],
    scorer: ScoringMatrix,
) -> tuple[np.ndarray, np.ndarray, list[tuple[int, tuple[int, ...]]]]:
    """Build N-dimensional score and direction grids.

    Args:
        seqs: Padded sequences (each with leading gap already prepended).
        scorer: Scoring matrix.

    Returns:
        Tuple of (score_grid, direction_grid, moves) where direction_grid
        stores bitmasks of co-optimal moves at each cell.
    """
    n = len(seqs)
    shape = tuple(len(s) for s in seqs)
    moves = make_moves(n)

    s_grid = np.full(shape, -np.inf, dtype=np.float64)
    d_grid = np.zeros(shape, dtype=np.uint64)

    # Origin
    origin = tuple(0 for _ in range(n))
    s_grid[origin] = 0.0

    for idx in np.ndindex(*shape):
        if idx == origin:
            continue

        best_score = -np.inf
        best_mask = 0

        for bitmask, delta in moves:
            # Compute predecessor
            pred = tuple(idx[d] + delta[d] for d in range(n))

            # Check validity: all coords >= 0
            if any(p < 0 for p in pred):
                continue

            # Build column symbols
            symbols = _column_symbols(seqs, idx, bitmask, scorer.gap)

            # Skip all-gap columns
            if all(s == scorer.gap for s in symbols):
                continue

            score = s_grid[pred] + scorer[symbols]

            if score > best_score + _EPS:
                best_score = score
                best_mask = bitmask
            elif abs(score - best_score) < _EPS:
                best_mask |= bitmask

        s_grid[idx] = best_score
        d_grid[idx] = best_mask

    return s_grid, d_grid, moves


def ndim_nw_backtrace(
    seqs: list[list[Hashable]],
    d_grid: np.ndarray,
    moves: list[tuple[int, tuple[int, ...]]],
    gap: Hashable,
    idx: tuple[int, ...] | None = None,
) -> list[list[list[Hashable]]]:
    """Backtrace through the N-dimensional direction grid.

    Args:
        seqs: Padded sequences (with leading gap).
        d_grid: Direction grid (bitmask of optimal moves per cell).
        moves: Move list from make_moves().
        gap: Gap symbol.
        idx: Starting cell (default: terminal cell).

    Returns:
        List of alignments, each a list of N sequences (lists of symbols).
    """
    n = len(seqs)
    if idx is None:
        idx = tuple(len(s) - 1 for s in seqs)

    origin = tuple(0 for _ in range(n))

    if idx == origin or any(i < 0 for i in idx):
        return [[[] for _ in range(n)]]

    # Collect which moves are set in the direction bitmask
    mask = int(d_grid[idx])
    active_moves = [(bitmask, delta) for bitmask, delta in moves if mask & bitmask]

    if not active_moves:
        return [[[] for _ in range(n)]]

    all_results = []
    for bitmask, delta in active_moves:
        pred = tuple(idx[d] + delta[d] for d in range(n))
        if any(p < 0 for p in pred):
            continue
        sub_results = ndim_nw_backtrace(seqs, d_grid, moves, gap, pred)
        for sub in sub_results:
            for d in range(n):
                if bitmask & (1 << d):
                    sub[d].append(seqs[d][idx[d]])
                else:
                    sub[d].append(gap)
            all_results.append(sub)

    return all_results


def ndim_nw_align(
    seqs: list[list[Hashable]],
    matrix: ScoringMatrix,
    k: int | None = None,
) -> list[Alignment]:
    """Perform true N-dimensional Needleman-Wunsch alignment.

    Args:
        seqs: Sequences to align (raw, without leading gap).
        matrix: Scoring matrix.
        k: Maximum number of alignments to return (default: all).

    Returns:
        Sorted list of best alignments, truncated to k.
    """
    # Prepend gap to each sequence (matching anw.py pattern)
    padded = [[matrix.gap, *list(s)] for s in seqs]

    _, d_grid, moves = ndim_nw_grids(padded, matrix)

    raw_alms = ndim_nw_backtrace(padded, d_grid, moves, matrix.gap)

    # Deduplicate and wrap as Alignment objects
    seen: set[tuple[tuple[Hashable, ...], ...]] = set()
    alms = []
    for raw in raw_alms:
        key = tuple(tuple(s) for s in raw)
        if key in seen:
            continue
        seen.add(key)
        alms.append(Alignment(raw, score_alignment(raw, matrix)))

    result = sort_alignments(alms)
    return result[:k] if k is not None else result
