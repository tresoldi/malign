"""
Module with functions for the k-best pairwise alignment.
"""

# Import Python standard libraries
from itertools import islice
from typing import Hashable, List, Optional, Sequence, Tuple

# Import 3rd party libraries
import networkx as nx

# Import other modules
from .alignment import Alignment
from .scoring_matrix import ScoringMatrix
from .utils import identity_matrix, pairwise_iter, score_alignment, sort_alignments


def compute_graph(
        seq_a: Sequence[Hashable], seq_b: Sequence[Hashable], matrix: ScoringMatrix
) -> nx.DiGraph:
    """
    Computes a weighted directed graph for alignment.

    The function builds a weighted directed graph for a pairwise alignment
    between two sequences, in a manner analogous to a Needleman-Wunsch
    alignment matrix. Alignment sites are used as graph nodes and
    transitions are used as weighted edges.

    Unlike in Needleman-Wunsch, no prior decision of transition "direction"
    is performed (i.e., no selection for the maximum in each cell),
    meaning that all possible diagonal, vertical, and horizontal directions
    are computed.

    As this alignment method is designed as a building block for the
    alignment of sequences with different alphabets and with non symmetric
    transitions costs, the alignment and the scorer are intended to align
    the second sequence (`seq_b`) in terms of the first one (`seq_a`).
    Sequences can include gap symbols.

    @param seq_a: A list of hashable elements to be aligned. The pairwise alignment
        is build in terms of this first sequence.
    @param seq_b: A list of hashable elements to be aligned.
    @param matrix: A ScoringMatrix the alignments. Dictionary keys
        are tuples in the format (`seq_a character`, `seq_b character`),
        dictionary values are numbers holding the
        score (the higher the score, the more favoured will the
        alignment be). Gaps have individual scores in relation to sequence
        values. The scorer itself is optional and, if not provided,
        a default one will be built. The implementation assumes the
        matrix has already been filled or will be automatically filled if
        necessary (with inference of missing values).
    @return: A Diretional Graph with alignment sites as nodes and transitions
        as weighted edges.
    """

    # Create a scorer if it was not provided and extract its maximum
    # value (score), so that scores are later corrected to costs when
    # adding edges.
    # Please note that, while it would make more sense or at least be more
    # natural, given our graph approach, to store *cost* values instead of
    # *score* values (also favoring smaller/negative ones), the "tradition" in
    # sequence alignment is to report scorer.
    max_score = max(matrix.scores.values())

    # Add gaps to the beginning of both sequences, emulating the first row
    # and first column in NW.
    # NOTE: the "+" operation on lists here allows us to, in a single step,
    # add the necessary alignment gap and make an in-memory copy of both
    # sequences, preserving the original ones.
    seq_a = [matrix.gap] + list(seq_a)
    seq_b = [matrix.gap] + list(seq_b)

    # Build the directional graph and add edges iterating from the bottom
    # right to the top left corner (as in NW).
    # The `networkx` library will automatically add the nodes. Note that, as intended
    # to our library, `symbol_a` and `symbol_b` might be gaps themselves.
    graph = nx.DiGraph()
    for i in range(len(seq_a) - 1, -1, -1):
        # Cache symbol
        symbol_a = seq_a[i]

        for j in range(len(seq_b) - 1, -1, -1):
            # Cache symbol
            symbol_b = seq_b[j]

            # Get the correspondence for the current cell (diagonal),
            # horizontal (gap in the top sequence, `seq_b`) and vertical
            # (gap in the left sequence, `seq_a`). We also check if we are
            # at a border, where diagonal and either vertical or horizontal
            # movement is not possible. Note that the logic is more
            # explicit to save computation time.
            if i == 0 and j == 0:
                dig_score = None
                hor_score = None
                ver_score = None
            elif i == 0:
                dig_score = None
                hor_score = None
                ver_score = matrix[symbol_a, matrix.gap]
            elif j == 0:
                dig_score = None
                hor_score = matrix[matrix.gap, symbol_b]
                ver_score = None
            else:
                dig_score = matrix[symbol_a, symbol_b]
                hor_score = matrix[matrix.gap, symbol_b]
                ver_score = matrix[symbol_a, matrix.gap]

            # Add edges (and nodes automatically); note that this, as expected, will
            # not add or adjust the top left cell (0, 0).
            if dig_score is not None:
                graph.add_edge((i - 1, j - 1), (i, j), weight=max_score - dig_score)
            if hor_score is not None:
                graph.add_edge((i - 1, j), (i, j), weight=max_score - hor_score)
            if ver_score is not None:
                graph.add_edge((i, j - 1), (i, j), weight=max_score - ver_score)

    return graph


def build_align(
        path: List[Tuple[int, int]],
        seq_a: Sequence[Hashable],
        seq_b: Sequence[Hashable],
        gap: Hashable = "-",
) -> Tuple[Sequence[Hashable], Sequence[Hashable]]:
    """
    Builds a pairwise alignment from a path of sequence indexes.

    The function receives a `path` as provided from a k shortest path
    searching method and applies it to a pair of sequences, building the
    alignment with the actual values and adding gaps if necessary.
    Paths do not need to start at the top left corner nor end at the bottom
    right one. Sequences are allowed to already have gaps in them.

    @param path: A list of string labels in the format `"seq_a_idx:seq_b_idx"`.
    @param seq_a: A list of elements used as the first sequence (on the vertical
        border).
    @param seq_b: A list of elements used as the second sequences (on the horizontal
        border).
    @param gap: A string used as gap symbol. Defaults to `"-"`.
    @return: The alignments for sequence A and B.
    """

    # We make a copy of the sequences so that we can pop from them without
    # consuming the original value. This also simplifies the code, as
    # we can slice from the original point in the path (allowing a starting
    # path different from (0,0).
    seq_a = list(seq_a[path[0][0]:])
    seq_b = list(seq_b[path[0][1]:])

    # Build the alignment sequences, adding (more) gaps if necessary
    alm_a = []
    alm_b = []
    for source, target in pairwise_iter(path):
        if target[1] == source[1]:
            # vertical movement
            alm_a.append(seq_a.pop(0))
            alm_b.append(gap)
        elif target[0] == source[0]:
            # horizontal movement
            alm_a.append(gap)
            alm_b.append(seq_b.pop(0))
        else:
            # diagonal movement
            alm_a.append(seq_a.pop(0))
            alm_b.append(seq_b.pop(0))

    return alm_a, alm_b


def align(
        graph: nx.DiGraph,
        ne_loc: Tuple[int, int],
        sw_loc: Tuple[int, int],
        seq_a: Sequence[Hashable],
        seq_b: Sequence[Hashable],
        matrix: ScoringMatrix,
        n_paths: Optional[int] = None,
) -> List[Alignment]:
    """
    Return the `k` best alignments in terms of costs.

    The function with take a path as source and target nodes in a graph,
    score the best alignments in terms of transition weights and gap
    penalties, and return the `k` best alignments. As an NP hard problem,
    it is not mathematically guaranteed that the best path will be found
    unless the search pass, defined by `n_paths`, includes all possible
    paths.

    @param graph: A directed weighted graph for the alignment, as returned by
        compute_graph().
    @param ne_loc:
    @param sw_loc:
    @param seq_a: A list of elements used as the first sequence (on the vertical
        border).
    @param seq_b: A list of elements used as the second sequences (on the horizontal
        border).
    @param matrix:
    @param n_paths: The number of alignment paths to be collected for scoring. If
        not provided, a default value will be calculated based on the
        length of the sequences and the complexity of the graph.
    @return: A sorted list of best alignments in the format `[alm, cost]`,
        with `alm` as a list of aligned sequences and `cost` the
        computed alignment cost including transitions and penalties.
    """

    # Given that `networkx` does not return the sum of edge weights and
    # that we need to perform individual score adjustments for
    # gap opening and gap extension, we first need to collect the `n`
    # shortest simple paths, build the alignments, and correct the scores
    # before yielding the top `k` alignments. The number of `n` paths
    # to collect is difficult to determine beforehand, so whether the
    # user is allowed to provide it or we determine it with this simple
    # operation.
    n_paths = n_paths or min(len(seq_a), len(seq_b))

    # Compute the paths from `networkx`; note that this is an iterator,
    # which is sliced and converted to a list in the loop below
    paths = nx.shortest_simple_paths(graph, ne_loc, sw_loc, weight="weight")

    # Iterate over the collected paths, collecting the corresponding
    # alignments and weights so we can sort after the loop
    alignments = []
    for path in list(islice(paths, n_paths)):
        # Build sequential representation of the alignment alignment
        alm_seq_a, alm_seq_b = build_align(path, seq_a, seq_b, gap=matrix.gap)
        score = score_alignment([alm_seq_a, alm_seq_b], matrix)
        alignments.append(Alignment([alm_seq_a, alm_seq_b], score))

    # Sort and return
    return sort_alignments(alignments)


def yenksp_align(
        seq_a: Sequence[Hashable],
        seq_b: Sequence[Hashable],
        k: Optional[int] = 4,
        matrix: Optional[ScoringMatrix] = None,
        ne_loc: Tuple[int, int] = (0, 0),
        sw_loc: Optional[Tuple[int, int]] = None,
) -> List[Alignment]:
    """
    Perform pairwise alignment with the Yen K Shortest Paths method.

    @param seq_a: The first sequence to be aligned.
    @param seq_b: The second sequence to be aligned.
    @param k: The number of paths to compute. The method can compute all possible paths
        by explicitly passing `None` as a value, but this might be prohibitively
        expansive depending on the length of the sequences and the complexity
        of the alphabets/matrix. Defaults to 4 (that is, twice the number of
        sequences).
    @param matrix: The matrix for the asymmetric scoring. Note that the order of the domains must
        follow the order of the sequences provided, that is, the matrix should be
        addressed with `(seq_a_symbol, seq_b_symbol)`. The implementation assumes the
        matrix has already been filled or will be automatically filled if
        necessary (with inference of missing values).
    @param ne_loc:
    @param sw_loc:
    @return:
    """

    # Obtain parameters
    if not matrix:
        matrix = identity_matrix([seq_a, seq_b])

    # Get the paths extremes, if not provided
    sw_loc = sw_loc or (len(seq_a), len(seq_b))

    # Compute the graph for the alignment
    graph = compute_graph(seq_a, seq_b, matrix)

    # Align from the graph and return
    return align(graph, ne_loc, sw_loc, seq_a, seq_b, matrix, n_paths=k)
