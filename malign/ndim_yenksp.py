"""N-dimensional Yen K-Shortest Paths alignment."""

from collections.abc import Hashable
from itertools import pairwise
from math import prod

import numpy as np

from .alignment import Alignment
from .ndim_common import _column_symbols, make_moves
from .scoring_matrix import ScoringMatrix
from .utils import pad_sequence, score_alignment, sort_alignments
from .yenksp import _add_edge, _new_graph, yen_ksp

# Maximum number of graph nodes before refusing to build.
# 10M nodes ≈ several GB of RAM for the adjacency structure.
MAX_GRAPH_NODES = 10_000_000


def compute_graph_nd(
    seqs: list[list[Hashable]],
    matrix: ScoringMatrix,
) -> dict:
    """Build N-dimensional weighted directed graph for alignment.

    Nodes are N-tuples of indices. Edge weights use the max_score - score
    conversion to turn maximization into shortest-path minimization.

    Args:
        seqs: Padded sequences (each with leading gap prepended).
        matrix: Scoring matrix.

    Returns:
        Graph as adjacency dict suitable for yen_ksp.

    Raises:
        MemoryError: If the grid would exceed MAX_GRAPH_NODES.
    """
    n = len(seqs)
    shape = tuple(len(s) for s in seqs)

    grid_size = prod(shape)
    if grid_size > MAX_GRAPH_NODES:
        dims_str = " × ".join(str(s) for s in shape)
        raise MemoryError(
            f"N-dimensional alignment grid too large: {dims_str} = {grid_size:,} nodes "
            f"(limit is {MAX_GRAPH_NODES:,}). Reduce the number or length of sequences."
        )

    moves = make_moves(n)
    max_score = max(matrix.scores.values())

    graph = _new_graph()

    for idx in np.ndindex(*shape):
        for bitmask, delta in moves:
            # Compute predecessor
            pred = tuple(idx[d] + delta[d] for d in range(n))
            if any(p < 0 for p in pred):
                continue

            symbols = _column_symbols(seqs, idx, bitmask, matrix.gap)
            if all(s == matrix.gap for s in symbols):
                continue

            score = matrix[symbols]
            weight = max_score - score
            _add_edge(graph, pred, idx, weight)

    return graph


def build_align_nd(
    path: list[tuple[int, ...]],
    seqs: list[list[Hashable]],
    gap: Hashable,
) -> list[list[Hashable]]:
    """Reconstruct N-dimensional alignment from a graph path.

    For each step, if a dimension's coordinate changed, the corresponding
    sequence symbol is used; otherwise a gap is inserted.

    Args:
        path: List of N-tuple nodes from source to target.
        seqs: Padded sequences (with leading gap).
        gap: Gap symbol.

    Returns:
        List of N aligned sequences.
    """
    n = len(seqs)
    # Track position in each sequence (skip the leading gap at index 0)
    pos = [path[0][d] for d in range(n)]

    result = [[] for _ in range(n)]
    for source, target in pairwise(path):
        for d in range(n):
            if target[d] != source[d]:
                result[d].append(seqs[d][target[d]])
                pos[d] = target[d]
            else:
                result[d].append(gap)

    return result


def ndim_yenksp_align(
    seqs: list[list[Hashable]],
    k: int,
    matrix: ScoringMatrix,
) -> list[Alignment]:
    """Perform true N-dimensional alignment using Yen's K-Shortest Paths.

    Args:
        seqs: Sequences to align (raw, without leading gap).
        k: Number of best alignments to find.
        matrix: Scoring matrix.

    Returns:
        Sorted list of best alignments.
    """
    padded = [pad_sequence(s, matrix.gap) for s in seqs]

    graph = compute_graph_nd(padded, matrix)

    n = len(seqs)
    origin = tuple(0 for _ in range(n))
    terminal = tuple(len(s) - 1 for s in padded)

    paths = yen_ksp(graph, origin, terminal, k)

    seen: set[tuple[tuple[Hashable, ...], ...]] = set()
    alms = []
    for path in paths:
        raw = build_align_nd(path, padded, matrix.gap)
        key = tuple(tuple(s) for s in raw)
        if key in seen:
            continue
        seen.add(key)
        alms.append(Alignment(raw, score_alignment(raw, matrix)))

    return sort_alignments(alms)
