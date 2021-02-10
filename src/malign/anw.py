"""
Module for computing asymmetric Needleman–Wunsch alignments.
"""

# Import Python standard libraries
import itertools
from typing import Optional, Hashable, Sequence, List, Tuple, Dict

# Import from package
from .utils import score_alignment, sort_alignments
from .scoring_matrix import ScoringMatrix

# Defines the map for directions; keys are tuples of three booleans,
# resulting from comparison of direction scores with the best scores,
# being, in order, (1) diagonal, (2) horizontal, (3) vertical. Used for
# debugging.
DIRECTION_MAP = {
    (False, False, False): "-",  # stationary, up left corner
    (True, True, True): "*",  # all movements with the same score
    (True, False, False): "↖",  # diagonal
    (False, True, False): "←",  # horizontal
    (False, False, True): "↑",  # vertical
    (True, True, False): "←↖",  # diagonal and horizontal
    (True, False, True): "↑↖",  # diagonal and vertical
    (False, True, True): "←↑",  # horizontal and vertical
}


def nw_grids(
    seq_a: List[Hashable], seq_b: List[Hashable], scorer: ScoringMatrix
) -> Tuple[List[List[float]], List[List[Tuple[bool, bool, bool]]]]:
    """
    Build the Needleman-Wunsch grids

    Note that the sequences must already have the initial gap added to them at this
    point.

    @param seq_a: First sequence for the pairwise alignment.
    @param seq_b: Second sequence for the pairwise alignment.
    @param scorer: The scoring matrix for the pairwise alignement.
    @return:
    """

    # cache lengths
    len_a, len_b = len(seq_a), len(seq_b)

    # Initialize (seq_a x seq_b) grids, one for the scores (`s_grid`) and one
    # for the directions (`d_grid`). Note that `seq_a` is modelled at the top
    # (i.e., columns), so the indexing is performed with `grid[b][a]`
    s_grid: List[List[float]] = [[0.0] * len_a for _ in seq_b]
    d_grid: List[List[Tuple[bool, bool, bool]]] = [
        [(False, False, False)] * len_a for _ in seq_b
    ]

    # Fill first row and column of both grids
    s_grid[0][0] = scorer[scorer.gap, scorer.gap]
    d_grid[0][0] = (False, False, False)  # no movement
    for i in range(1, len_a):
        s_grid[0][i] = -i
        d_grid[0][i] = (False, True, False)  # horizontal
    for j in range(1, len_b):
        s_grid[j][0] = -j
        d_grid[j][0] = (False, False, True)  # vertical

    # Fill the cells by (a) computing diagonal, vertical, and horizontal cost,
    # (b) getting the highest value for `s_grid`, (c) checking the directions
    # that match such highest value for `d_grid`
    for i, j in itertools.product(range(1, len_a), range(1, len_b)):
        # compute direction scorers
        diag = s_grid[j - 1][i - 1] + scorer[seq_a[i], seq_b[j]]
        horz = s_grid[j][i - 1] + scorer[seq_a[i], scorer.gap]
        vert = s_grid[j - 1][i] + scorer[scorer.gap, seq_b[j]]

        # get best score and matching tuple
        best_score = max([diag, horz, vert])
        match_dir = (diag == best_score, horz == best_score, vert == best_score)

        # set values
        s_grid[j][i] = best_score
        d_grid[j][i] = match_dir

    return s_grid, d_grid


def _nw_product(
    prev_alms: List[Dict[str, List[Hashable]]],
    char_a: Hashable,
    char_b: Hashable,
    paths: List[Dict[Hashable, List[Hashable]]],
) -> List[Dict[str, List[Hashable]]]:
    """
    Internal function for building a product of paths.

    The function is used for NW alignments with two or more directions.
    """
    ret_alms = []
    for alm in prev_alms:
        ret_alms += [
            {"a": [*alm["a"], char_a, *path["a"]], "b": [*alm["b"], char_b, *path["b"]]}
            for path in paths
        ]

    return ret_alms


# pylint: disable=too-many-branches
def nw_backtrace(
    seq_a: List[Hashable],
    seq_b: List[Hashable],
    d_grid: List[List[Tuple[bool, bool, bool]]],
    gap: Hashable,
    i: Optional[int] = None,
    j: Optional[int] = None,
) -> List[Dict[Hashable, List[Hashable]]]:
    """
    Run the Needleman-Wunsch backtrace operation.

    Note that the alignments are returned in reverse order, from the bottom right to the
    top left; as the function is called recursively, it is up to the function
    caller to reverse it (if so desired).

    As this function will only operate pairwise, for easiness of debugging and
    inspection the returned structure does not follow the approach used elsewhere of
    a list of sequences, but it is a dictionary with `a` and `b` keys.

    @param seq_a: First sequence for the pairwise alignment.
    @param seq_b: Second sequence for the pairwise alignment.
    @param d_grid: A direction grid, as returned from `nw_grids()`.
    @param gap:
    @param i:
    @param j:
    @return:
    """

    # Get parameters, and default for the full alignment, if (i, j) is not provided
    if not i and not j:
        i = len(seq_a) - 1
        j = len(seq_b) - 1

    # Define empty, initial alignment collection with a single alignment
    # TODO: as it is pairwise, we don´t really need "a" and "b" and can index a normal list
    alms = [{"a": [], "b": []}]

    # Does the backtrace, using recursion when necessary and changing in
    # place as much as possible for speed
    while True:
        if d_grid[j][i] == (True, False, False):
            # diagonal
            for alm in alms:
                alm["a"].append(seq_a[i])
                alm["b"].append(seq_b[j])
            i, j = i - 1, j - 1
        elif d_grid[j][i] == (False, True, False):
            # horizontal
            for alm in alms:
                alm["a"].append(seq_a[i])
                alm["b"].append(gap)
            i = i - 1
        elif d_grid[j][i] == (False, False, True):
            # vertical
            for alm in alms:
                alm["a"].append(gap)
                alm["b"].append(seq_b[j])
            j = j - 1
        elif d_grid[j][i] == (True, False, True):
            # diagonal and vertical
            diag_paths = nw_backtrace(seq_a, seq_b, d_grid, gap, i - 1, j - 1)
            vert_paths = nw_backtrace(seq_a, seq_b, d_grid, gap, i, j - 1)

            ret_alms = _nw_product(alms, seq_a[i], seq_b[j], diag_paths)
            ret_alms += _nw_product(alms, gap, seq_b[j], vert_paths)

            return ret_alms
        elif d_grid[j][i] == (True, True, False):
            # diagonal and horizontal
            diag_paths = nw_backtrace(seq_a, seq_b, d_grid, gap, i - 1, j - 1)
            horz_paths = nw_backtrace(seq_a, seq_b, d_grid, gap, i - 1, j)

            ret_alms = _nw_product(alms, seq_a[i], seq_b[j], diag_paths)
            ret_alms += _nw_product(alms, seq_a[i], gap, horz_paths)

            return ret_alms
        elif d_grid[j][i] == (False, True, True):
            # vertical and horizontal
            vert_paths = nw_backtrace(seq_a, seq_b, d_grid, gap, i, j - 1)
            horz_paths = nw_backtrace(seq_a, seq_b, d_grid, gap, i - 1, j)

            ret_alms = _nw_product(alms, gap, seq_b[j], vert_paths)
            ret_alms += _nw_product(alms, seq_a[i], gap, horz_paths)

            return ret_alms
        elif d_grid[j][i] == (True, True, True):
            # diagonal, vertical, and horizontal
            diag_paths = nw_backtrace(seq_a, seq_b, d_grid, gap, i - 1, j - 1)
            vert_paths = nw_backtrace(seq_a, seq_b, d_grid, gap, i, j - 1)
            horz_paths = nw_backtrace(seq_a, seq_b, d_grid, gap, i - 1, j)

            ret_alms = _nw_product(alms, seq_a[i], seq_b[j], diag_paths)
            ret_alms += _nw_product(alms, gap, seq_b[j], vert_paths)
            ret_alms += _nw_product(alms, seq_a[i], gap, horz_paths)

            return ret_alms
        else:
            # NOTE: this exception would never be reached in normal and "correct"
            # operation, but it is needed here as the user manipulation of the grids
            # and/or matrices, which is intentional for advanced cases, could mean
            # stationary cells
            raise ValueError(f"Missing direction {d_grid[j][i]} at (i={i}, j={j})")

        # Break when reaching the top left corner
        if i == 0 and j == 0:
            break

    return alms


def nw_align(
    seq_a: Sequence[Hashable],
    seq_b: Sequence[Hashable],
    matrix: ScoringMatrix,
    k: Optional[int] = None,
):
    """
    Perform pairwise alignment with the Asymmetric Needleman-Wunsch method.

    @param seq_a: The first sequence to be aligned.
    @param seq_b: The second sequence to be aligned.
    @param matrix: The matrix for the asymmetric scoring. Note that the order of the domains must
        follow the order of the sequences provided, that is, the matrix should be
        addressed with `(seq_a_symbol, seq_b_symbol)`. The implementation assumes the
        matrix has already been filled or will be automatically filled if
        necessary (with inference of missing values).
    @param k: Number of alignments to include in return. Note that is the upper limit, as it
        impossible to guarantee that there will be as many alignments as requested (due
        to both the sequences and the scoring matrix) and that this implementation of
        the Needleman-Wunsch algorithm is not intended for collection of as many
        k-best alignments as possible (if so desired, the `yenksp` method is
        recommended). Defaults to `None`, meaning that all the collected alignments
        will be returned.
    @return:
    """

    # Add initial gaps; note that this also makes a
    # copy of the contents of each sequence, so we preserve the original
    # memory in-place
    seq_a = [matrix.gap] + list(seq_a)
    seq_b = [matrix.gap] + list(seq_b)

    # Build Needleman-Wunsch grids; note that the scoring grid (the first value returned
    # by `nw_grids()`) is not used in this routine, as the scoring is performed with the
    # more complete `score_alignment()` function.
    # TODO: do we really need to pass matrix.gap instead of reading it from there? the method is more general here...
    _, d_grid = nw_grids(seq_a, seq_b, matrix)

    # Obtain the alignments from backtrace, and them along with a score;
    # sequences are reversed after it with `[::-1]`, as NW follows a northwest-direction;
    # this is not necessary for computing the score, as it is performed per site and
    # thus the mirrored sequences have no effect
    alms = [
        {
            "seqs": [alm["a"][::-1], alm["b"][::-1]],
            "score": score_alignment([alm["a"], alm["b"]], matrix),
        }
        for alm in nw_backtrace(seq_a, seq_b, d_grid, matrix.gap)
    ]

    # Sort and return
    return sort_alignments(alms)[:k]
