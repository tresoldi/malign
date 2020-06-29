"""
Module for computing Needleman–Wunsch alignments.
"""

# Import Python standard libraries
import itertools

# Import third-party libraries
import numpy as np

# Import other modules
import malign.utils as utils

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

# NOTE seq_a and seq_b must already have the initial gap
# TODO: pass gap symbol
def nw_grids(seq_a, seq_b, scorer):
    gap = "-"

    # cache lengths
    len_a = len(seq_a)
    len_b = len(seq_b)

    # Initialize (seq_a x seq_b) grids, one for the scores (`s_grid`) and one
    # for the directions (`d_grid`). Note that `seq_a` is modelled at the top
    # (i.e., columns), so the indexing is performed with `grid[b][a]`
    s_grid = np.empty([len_b, len_a])
    d_grid = np.empty([len_b, len_a])

    # Fill first row and column of both grids
    s_grid[0][0] = scorer["-", "-"]
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


# TODO: implement blocks
# TODO: allows collections of i_stop, j_stop?
# TODO: default starting i/j to the bottom right
# TODO: add gap
# TODO: note that the alignments will be returned in reverse order; as
#       the function is called recursively, it is up to the fuctional call
#       ing it to reverse it back (if desired)
# NOTE: allowing all combinations of directions, even those theoretically
#       impossible with NW, because we will allow for manual changes
def nw_backtrace(seq_a, seq_b, d_grid, i, j, i_stop=0, j_stop=0):
    gap = "-"

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
            # TODO: make sure code never gets stationary and remove this
            # TODO: note that loop could actually reach this if manually set
            #       -- deal by just returning
            err = "missing_dir %i at (i=%i, j=%i)" % (d_grid[j][i], i, j)
            raise ValueError(err)

        if i == i_stop and j == j_stop:
            break

    return alms


# Score an alignment with standard NW method (symmetric, etc)
# TODO: support for different scoring for gaps at extremeties
# TODO: gap symbol, gap opening, gap penalty
# TODO: rename to nw_add_scores to keep consisntency in namespace
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
        # TODO: mean of score_a and score_b?
        alm["score_a"] = (sum(gaps_a) * gap_pen) + (len(gaps_a) * gap_opn)
        alm["score_b"] = (sum(gaps_b) * gap_pen) + (len(gaps_b) * gap_opn)
        alm["score"] = alm["score_a"] + alm["score_b"]


# TODO: pass gap
def nw_align(seq_a, seq_b, gap="-", scorer=None, k=1):
    # Build normal scorer if not provided
    if not scorer:
        scorer = utils.build_basic_scorer(set(seq_a), set(seq_b))

    # Add initial gaps; note that this also makes a
    # copy of the contents of each sequence, so we preserve the original
    # memory in-place
    seq_a = [gap] + seq_a
    seq_b = [gap] + seq_b

    # Build nw grids
    s_grid, d_grid = nw_grids(seq_a, seq_b, scorer)

    # move in the direction grid, reversing sequences after it
    i, j = len(seq_a) - 1, len(seq_b) - 1
    alms = nw_backtrace(seq_a, seq_b, d_grid, i, j)
    alms = [{"a": alm["a"][::-1], "b": alm["b"][::-1]} for alm in alms]
    add_nw_scores(alms)

    # Sort and return
    # TODO: checks, test, better order, etc.
    # TODO: make general utils function
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
