"""Module for computing asymmetric Needleman-Wunsch alignments."""

import itertools
from collections.abc import Hashable, Sequence

from .alignment import Alignment
from .scoring_matrix import ScoringMatrix
from .utils import score_alignment, sort_alignments

# Direction map: keys are (diagonal, horizontal, vertical) match tuples.
DIRECTION_MAP = {
    (False, False, False): "-",
    (True, True, True): "*",
    (True, False, False): "↖",
    (False, True, False): "←",
    (False, False, True): "↑",
    (True, True, False): "←↖",
    (True, False, True): "↑↖",
    (False, True, True): "←↑",
}


def nw_grids(
    seq_a: list[Hashable],
    seq_b: list[Hashable],
    scorer: ScoringMatrix,
) -> tuple[list[list[float]], list[list[tuple[bool, bool, bool]]]]:
    """Build the Needleman-Wunsch score and direction grids.

    Sequences must already have the initial gap prepended.

    Args:
        seq_a: First sequence (with leading gap).
        seq_b: Second sequence (with leading gap).
        scorer: Scoring matrix for the pairwise alignment.

    Returns:
        Tuple of (score_grid, direction_grid).
    """
    len_a, len_b = len(seq_a), len(seq_b)

    s_grid: list[list[float]] = [[0.0] * len_a for _ in seq_b]
    d_grid: list[list[tuple[bool, bool, bool]]] = [[(False, False, False)] * len_a for _ in seq_b]

    s_grid[0][0] = scorer[scorer.gap, scorer.gap]
    d_grid[0][0] = (False, False, False)
    for i in range(1, len_a):
        s_grid[0][i] = -i
        d_grid[0][i] = (False, True, False)
    for j in range(1, len_b):
        s_grid[j][0] = -j
        d_grid[j][0] = (False, False, True)

    for i, j in itertools.product(range(1, len_a), range(1, len_b)):
        diag = s_grid[j - 1][i - 1] + scorer[seq_a[i], seq_b[j]]
        horz = s_grid[j][i - 1] + scorer[seq_a[i], scorer.gap]
        vert = s_grid[j - 1][i] + scorer[scorer.gap, seq_b[j]]

        best_score = max([diag, horz, vert])
        match_dir = (diag == best_score, horz == best_score, vert == best_score)

        s_grid[j][i] = best_score
        d_grid[j][i] = match_dir

    return s_grid, d_grid


def _nw_product(
    prev_alms: list[dict[str, list[Hashable]]],
    char_a: Hashable,
    char_b: Hashable,
    paths: list[dict[Hashable, list[Hashable]]],
) -> list[dict[str, list[Hashable]]]:
    """Build a product of paths for NW alignments with multiple directions."""
    ret_alms = []
    for alm in prev_alms:
        ret_alms += [
            {"a": [*alm["a"], char_a, *path["a"]], "b": [*alm["b"], char_b, *path["b"]]}
            for path in paths
        ]

    return ret_alms


def nw_backtrace(
    seq_a: list[Hashable],
    seq_b: list[Hashable],
    d_grid: list[list[tuple[bool, bool, bool]]],
    gap: Hashable,
    i: int | None = None,
    j: int | None = None,
) -> list[dict[Hashable, list[Hashable]]]:
    """Run the Needleman-Wunsch backtrace operation.

    Alignments are returned in reverse order (bottom-right to top-left);
    the caller is responsible for reversing them.

    Args:
        seq_a: First sequence (with leading gap).
        seq_b: Second sequence (with leading gap).
        d_grid: Direction grid from nw_grids().
        gap: Gap symbol.
        i: Starting column index (default: last column).
        j: Starting row index (default: last row).

    Returns:
        List of alignment dictionaries with 'a' and 'b' keys.
    """
    alms = [{"a": [], "b": []}]

    if not i and not j:
        i = len(seq_a) - 1
        j = len(seq_b) - 1

    while True:
        if d_grid[j][i] == (True, False, False):
            for alm in alms:
                alm["a"].append(seq_a[i])
                alm["b"].append(seq_b[j])
            i, j = i - 1, j - 1
        elif d_grid[j][i] == (False, True, False):
            for alm in alms:
                alm["a"].append(seq_a[i])
                alm["b"].append(gap)
            i = i - 1
        elif d_grid[j][i] == (False, False, True):
            for alm in alms:
                alm["a"].append(gap)
                alm["b"].append(seq_b[j])
            j = j - 1
        elif d_grid[j][i] == (True, False, True):
            diag_paths = nw_backtrace(seq_a, seq_b, d_grid, gap, i - 1, j - 1)
            vert_paths = nw_backtrace(seq_a, seq_b, d_grid, gap, i, j - 1)

            ret_alms = _nw_product(alms, seq_a[i], seq_b[j], diag_paths)
            ret_alms += _nw_product(alms, gap, seq_b[j], vert_paths)

            return ret_alms
        elif d_grid[j][i] == (True, True, False):
            diag_paths = nw_backtrace(seq_a, seq_b, d_grid, gap, i - 1, j - 1)
            horz_paths = nw_backtrace(seq_a, seq_b, d_grid, gap, i - 1, j)

            ret_alms = _nw_product(alms, seq_a[i], seq_b[j], diag_paths)
            ret_alms += _nw_product(alms, seq_a[i], gap, horz_paths)

            return ret_alms
        elif d_grid[j][i] == (False, True, True):
            vert_paths = nw_backtrace(seq_a, seq_b, d_grid, gap, i, j - 1)
            horz_paths = nw_backtrace(seq_a, seq_b, d_grid, gap, i - 1, j)

            ret_alms = _nw_product(alms, gap, seq_b[j], vert_paths)
            ret_alms += _nw_product(alms, seq_a[i], gap, horz_paths)

            return ret_alms
        elif d_grid[j][i] == (True, True, True):
            diag_paths = nw_backtrace(seq_a, seq_b, d_grid, gap, i - 1, j - 1)
            vert_paths = nw_backtrace(seq_a, seq_b, d_grid, gap, i, j - 1)
            horz_paths = nw_backtrace(seq_a, seq_b, d_grid, gap, i - 1, j)

            ret_alms = _nw_product(alms, seq_a[i], seq_b[j], diag_paths)
            ret_alms += _nw_product(alms, gap, seq_b[j], vert_paths)
            ret_alms += _nw_product(alms, seq_a[i], gap, horz_paths)

            return ret_alms
        else:
            raise ValueError(f"Missing direction {d_grid[j][i]} at (i={i}, j={j})")

        if i == 0 and j == 0:
            break

    return alms


def nw_align(
    seq_a: Sequence[Hashable],
    seq_b: Sequence[Hashable],
    matrix: ScoringMatrix,
    k: int | None = None,
) -> list[Alignment]:
    """Perform pairwise alignment with the Asymmetric Needleman-Wunsch method.

    Args:
        seq_a: First sequence to align.
        seq_b: Second sequence to align.
        matrix: Scoring matrix (addressed as (seq_a_symbol, seq_b_symbol)).
        k: Maximum number of alignments to return (default: all).

    Returns:
        Sorted list of best alignments.
    """
    seq_a = [matrix.gap, *list(seq_a)]
    seq_b = [matrix.gap, *list(seq_b)]

    _, d_grid = nw_grids(seq_a, seq_b, matrix)

    alms = [
        Alignment(
            [alm["a"][::-1], alm["b"][::-1]],
            score_alignment([alm["a"], alm["b"]], matrix),
        )
        for alm in nw_backtrace(seq_a, seq_b, d_grid, matrix.gap)
    ]

    return sort_alignments(alms)[:k]
