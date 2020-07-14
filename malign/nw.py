"""
Module for computing Needleman–Wunsch alignments.
"""

# Import Python standard libraries
import itertools

# Import third-party libraries
import numpy as np

# Defines the map for directions; keys are tuples of three booleans,
# resulting from comparison of direction scores with the best scores,
# being, in order, (1) diagonal, (2) horizontal, (3) vertical
DIRECTION_MAP = {
    (False, False, False): 0,  # stationary, up left corner
    (True, True, True): 1,  # all movements with the same score
    (True, False, False): 2,  # diagonal
    (False, True, False): 3,  # horizontal
    (False, False, True): 4,  # vertical
    (True, True, False): 5,  # diagonal and horizontal
    (True, False, True): 6,  # diagonal and vertical
    (False, True, True): 7,  # horizonta and vertical
}


def nw_grids(seq_a, seq_b, scorer, gap):
    """
    Build the Needleman-Wunsch grids

    Note that the sequences must already have the initial gap added to them at this
    point.
    """

    # cache lengths
    len_a, len_b = len(seq_a), len(seq_b)

    # Initialize (seq_a x seq_b) grids, one for the scores (`s_grid`) and one
    # for the directions (`d_grid`). Note that `seq_a` is modelled at the top
    # (i.e., columns), so the indexing is performed with `grid[b][a]`
    s_grid = np.empty([len_b, len_a])
    d_grid = np.empty([len_b, len_a])

    # Fill first row and column of both grids
    s_grid[0][0] = scorer[gap, gap]
    d_grid[0][0] = DIRECTION_MAP[False, False, False]  # no movement
    for i in range(1, len_a):
        s_grid[0][i] = -i
        d_grid[0][i] = DIRECTION_MAP[False, True, False]  # horizontal
    for j in range(1, len_b):
        s_grid[j][0] = -j
        d_grid[j][0] = DIRECTION_MAP[False, False, True]  # vertical

    # Fill the cells by (a) computing diagonal, vertical, and horizontal cost,
    # (b) getting the highest value for `s_grid`, (c) checking the directions
    # that match such highest value for `d_grid`
    for i, j in itertools.product(range(1, len_a), range(1, len_b)):
        # compute direction scorers
        diag = s_grid[j - 1][i - 1] + scorer[seq_a[i], seq_b[j]]
        horz = s_grid[j][i - 1] - 1
        vert = s_grid[j - 1][i] - 1

        # get best score and matching tuple
        best_score = max([diag, horz, vert])
        match_dir = [score == best_score for score in [diag, horz, vert]]

        # set values
        s_grid[j][i] = best_score
        d_grid[j][i] = DIRECTION_MAP[tuple(match_dir)]

    return s_grid, d_grid


# pylint: disable=too-many-branches
def nw_backtrace(seq_a, seq_b, d_grid, i=None, j=None, **kwargs):
    """
    Run the Needleman-Wunsch backtrace operation.

    Note that the alignments are returned in reverse order, from the bottom right to the
    top left; as the function is called recursively, it is up to the function
    caller to reverse it (if so desired).
    """

    # Get parameters, and default to the full alignment, if (i, j) is not provided
    gap = kwargs.get("gap", "-")
    if not i and not j:
        i = len(seq_a) - 1
        j = len(seq_b) - 1

    # Define empty, initial alignment collection with a single alignment
    alms = [{"a": [], "b": []}]

    # Build product of paths, for cases with two or more directions
    def _product(prev_alms, char_a, char_b, paths):
        ret_alms = []
        for alm in prev_alms:
            ret_alms += [
                {
                    "a": [*alm["a"], char_a, *path["a"]],
                    "b": [*alm["b"], char_b, *path["b"]],
                }
                for path in paths
            ]
        return ret_alms

    # Does the backtrace, using recursion when necessary and changing in
    # place as much as possible for speed
    while True:
        if d_grid[j][i] == DIRECTION_MAP[True, False, False]:
            # diagonal
            for alm in alms:
                alm["a"].append(seq_a[i])
                alm["b"].append(seq_b[j])
            i, j = i - 1, j - 1
        elif d_grid[j][i] == DIRECTION_MAP[False, True, False]:
            # horizontal
            for alm in alms:
                alm["a"].append(seq_a[i])
                alm["b"].append(gap)
            i = i - 1
        elif d_grid[j][i] == DIRECTION_MAP[False, False, True]:
            # vertical
            for alm in alms:
                alm["a"].append(gap)
                alm["b"].append(seq_b[j])
            j = j - 1
        elif d_grid[j][i] == DIRECTION_MAP[True, False, True]:
            # diagonal and vertical
            diag_paths = nw_backtrace(seq_a, seq_b, d_grid, i - 1, j - 1)
            vert_paths = nw_backtrace(seq_a, seq_b, d_grid, i, j - 1)

            ret_alms = _product(alms, seq_a[i], seq_b[j], diag_paths)
            ret_alms += _product(alms, gap, seq_b[j], vert_paths)

            return ret_alms
        elif d_grid[j][i] == DIRECTION_MAP[True, True, False]:
            # diagonal and horizontal
            diag_paths = nw_backtrace(seq_a, seq_b, d_grid, i - 1, j - 1)
            horz_paths = nw_backtrace(seq_a, seq_b, d_grid, i - 1, j)

            ret_alms = _product(alms, seq_a[i], seq_b[j], diag_paths)
            ret_alms += _product(alms, seq_a[i], gap, horz_paths)

            return ret_alms
        elif d_grid[j][i] == DIRECTION_MAP[False, True, True]:
            # vertical and horizontal
            vert_paths = nw_backtrace(seq_a, seq_b, d_grid, i, j - 1)
            horz_paths = nw_backtrace(seq_a, seq_b, d_grid, i - 1, j)

            ret_alms = _product(alms, gap, seq_b[j], vert_paths)
            ret_alms += _product(alms, seq_a[i], gap, horz_paths)

            return ret_alms
        elif d_grid[j][i] == DIRECTION_MAP[True, True, True]:
            # diagonal, vertical, and horizontal
            diag_paths = nw_backtrace(seq_a, seq_b, d_grid, i - 1, j - 1)
            vert_paths = nw_backtrace(seq_a, seq_b, d_grid, i, j - 1)
            horz_paths = nw_backtrace(seq_a, seq_b, d_grid, i - 1, j)

            ret_alms = _product(alms, seq_a[i], seq_b[j], diag_paths)
            ret_alms += _product(alms, gap, seq_b[j], vert_paths)
            ret_alms += _product(alms, seq_a[i], gap, horz_paths)

            return ret_alms
        else:
            # NOTE: this exception would never be reached in normal and "correct"
            # operation, but it is needed here as the user manipulation of the grids
            # and/or matrices, which is intentional for advanced cases, could mean
            # stationary cells
            raise ValueError(f"Missing direction {d_grid[j][i]} at (i={i}, j={j})")

        if i == 0 and j == 0:
            break

    return alms


# Score an alignment with standard NW method (symmetric, etc)
# TODO: support for different scoring for gaps at extremeties
# TODO: gap symbol, gap opening, gap penalty
def add_nw_scores(alms):
    gap = "-"
    gap_pen = -1
    gap_opn = -2

    # Add score to each alignment in `alms`
    for alm in alms:
        # Get the gap groups
        gaps_a = [
            len(list(group))
            for char, group in itertools.groupby(alm["a"])
            if char == gap
        ]
        gaps_b = [
            len(list(group))
            for char, group in itertools.groupby(alm["b"])
            if char == gap
        ]

        # Compute scores with gap opening and penalty
        alm["score_a"] = (sum(gaps_a) * gap_pen) + (len(gaps_a) * gap_opn)
        alm["score_b"] = (sum(gaps_b) * gap_pen) + (len(gaps_b) * gap_opn)
        alm["score"] = (alm["score_a"] + alm["score_b"]) / 2.0


def nw_align(seq_a, seq_b, matrix, gap="-", k=1):
    """
    Perform pairwise alignment with the Asymmetric Needleman-Wunsch method.

    Parameters
    ==========
    seq_a : list
        The first sequence to be aligned.
    seq_b : list
        The second sequence to be aligned.
    matrix : dict or ScoringMatrix
        The matrix for the asymmetric scoring. Note that the order of the domains must
        follow the order of the sequences provided, that is, the matrix should be
        addressed with [seq_a_symbol, seq_b_symbol]. The implementation assumes the
        matrix has already been filled or will be automatically filled if
        necessary (with inference of missing values).
    gap : string
        The symbol to be used for gap representation. If `matrix` is a ScoringMatrix,
        it must match the gap symbol specified in it. Defaults to `"-"`.
    k : int
        Number of alignments to include in return. Note that is the upper limit, as it
        impossible to guarantee that there will be as many alignments as requested (due
        to both the sequences and the scoring matrix) and that this implementation of
        the Needleman-Wunsch algorithm is not intended for collection of as many
        k-best alignments as possible (if so desired, the `yenksp` method is
        recommended). Defaults to `1`; if equal to `None`, all the collected alignments
        will be returned.
    """

    # Validate parameters
    if not gap:
        raise ValueError("Gap symbol must be a non-empty string.")
    if not isinstance(matrix, dict):  # ScoringMatrix
        if gap != matrix.gap:
            raise ValueError("Different gap symbols.")

    # Add initial gaps; note that this also makes a
    # copy of the contents of each sequence, so we preserve the original
    # memory in-place
    seq_a = [gap] + seq_a
    seq_b = [gap] + seq_b

    # Build Needleman-Wunsch grids; note that the scoring grid is not used in this
    # implementation.
    # pylint: disable=unused-variable
    s_grid, d_grid = nw_grids(seq_a, seq_b, matrix, gap)

    # move in the direction grid, reversing sequences after it with `[::-1]`
    alms = nw_backtrace(seq_a, seq_b, d_grid)
    alms = [{"a": alm["a"][::-1], "b": alm["b"][::-1]} for alm in alms]

    # TODO: use matrix
    add_nw_scores(alms)

    # Sort and return
    # TODO: checks, test, better order, etc.
    alms = sorted(
        alms,
        key=lambda alm: (
            -alm["score"],
            -alm["score_a"],
            -alm["score_b"],
            alm["a"],
            alm["b"],
        ),
    )

    return alms[:k]
