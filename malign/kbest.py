"""
Module with functions for the k-best pairwise alignment.
"""

# Import Python standard libraries
from itertools import product, islice, tee, groupby
import operator

# Import 3rd party libraries
import networkx as nx

# Import other modules
import malign.utils as utils

# Define an internal auxiliary function.
def _pairwise_iter(iterable):
    """
    Internal function for sequential pairwise iteration.
    The function follows the recipe in Python's itertools documentation.
    [https://docs.python.org/3/library/itertools.html]
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    """

    item_a, item_b = tee(iterable)
    next(item_a, None)

    return zip(item_a, item_b)


# TODO: allow different gap symbols
# TODO: check what happens with the ("-", "-") default of 0.0 when adjusting
#       weight for cost
# TODO: should add notes in SE direction, instead of NW? This could be
#       potentially be used to improve path finding, combining/selecting the
#       results of both directions, but it *cannot* be mixed in a single
#       graph, as we can't have loops for Yen's -- maybe this also implies
#       adding gaps to the end, in a mirrored needleman-wunsch
# TODO: function to export/represent the graph, or the matrix (should it be
#       an object?)
# TODO: allow to correct costs by some parameter, perhaps even non-linear?
# TODO: investigate usage of tuples of integers as keys
def compute_graph(seq_a, seq_b, matrix=None):
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

    Parameters
    ==========

    seq_a : list or iterable
        A list of hasheable elements to be aligned. The pairwise alignment
        is build in terms of this first sequence.
    seq_b : list or iterable
        A list of hasheable elements to be aligned.
    matrix: dict
        A dictionary of scoring weights for the alignments. Dictionary keys
        are tuples in the format (`seq_a character`, `seq_b character`),
        dictionary values are numbers holding the
        score (the higher the score, the more favoured will the
        alignment be). Gaps have individual scores in relation to sequence
        values. The scorer itself is optional and, if not provided,
        a default one will be built.

    Returns
    =======
    graph : networkx graph
        A Diretional Graph with alignment sites as nodes and transitions
        as weighted edges.
    """

    # Create a scorer if it was not provided and extract its maximum
    # value (score), so that scores are later corrected to costs when
    # adding edges.
    # Please note that, while it would make more sense or at least be more
    # natural, given our graph approach, to store *cost* values instead of
    # *score* vaues (also favoring smaller/negative ones), the "tradition" in
    # sequence alignment is to report scorer.
    # TODO: assuming matrix is filled
    #    scorer = utils.fill_matrix(set(seq_a), set(seq_b), scorer)
    max_score = max(matrix.scores.values())

    # Add gaps to the beginning of both sequences, emulating the first row
    # and first column in NW.
    # NOTE: the "+" operation on lists here allows us to, in a single step,
    # add the necessary alignment gap and make an in-memory copy of both
    # sequences, preserving the original ones.
    seq_a = ["-"] + list(seq_a)
    seq_b = ["-"] + list(seq_b)

    # Build the directional graph and add edges iterating from the bottom
    # right to the top left corner (as in NW).
    # NOTE: the `networkx` library will automatically add the nodes.
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
                ver_score = matrix[symbol_a, "-"]
            elif j == 0:
                dig_score = None
                hor_score = matrix["-", symbol_b]
                ver_score = None
            else:
                dig_score = matrix[symbol_a, symbol_b]
                hor_score = matrix["-", symbol_b]
                ver_score = matrix[symbol_a, "-"]

            # Add edges (and nodes automatically)
            if dig_score is not None:
                graph.add_edge(
                    "%i:%i" % (i - 1, j - 1),
                    "%i:%i" % (i, j),
                    weight=max_score - dig_score,
                )
            if hor_score is not None:
                graph.add_edge(
                    "%i:%i" % (i - 1, j), "%i:%i" % (i, j), weight=max_score - hor_score
                )
            if ver_score is not None:
                graph.add_edge(
                    "%i:%i" % (i, j - 1), "%i:%i" % (i, j), weight=max_score - ver_score
                )

    return graph


def build_align(path, seq_a, seq_b, gap="-"):
    """
    Builds a pairwise alignment from a path of sequence indexes.
    The function receives a `path` as provided from a k shortest path
    searching method and applies it to a pair of sequences, building the
    alignment with the actual values and adding gaps if necessary.
    Paths do not need to start at the top left corner nor end at the bottom
    right one. Sequences are allowed to already have gaps in them.

    Parameters
    ==========
    path : list
        A list of string labels in the format `"seq_a_idx:seq_b_idx"`.
    seq_a : list
        A list of elements used as the first sequence (on the vertical
        border).
    seq_b : list
        A list of elements used as the second sequences (on the horizontal
        border).
    gap : str
        A string used as gap symbol. Defaults to `"-"`.

    Returns
    =======
    alm_a, alm_b : lists
        The alignments for sequence A and B.
    """

    # Transform `path`, which is a list of string labels in the format
    # `"seq_a_idx:seq_b_idx"`, into an easily iterable list of tuples
    # in the format (seq_a_idx, seq_b_idx). We can then extract the
    # initial values of `prev_i` and `prev_j`.
    path = [[int(v) for v in cell.split(":")] for cell in path]

    # We make a copy of the sequences so that we can pop from them without
    # consuming the original value. This also simplifies the code, as
    # we can slice from the original point in the path (allowing a starting
    # path different from (0,0).
    seq_a = list(seq_a[path[0][0] :])
    seq_b = list(seq_b[path[0][1] :])

    # Build the alignment sequences, adding (more) gaps if necessary
    alm_a = []
    alm_b = []
    for source, target in _pairwise_iter(path):
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

    return {"a": alm_a, "b": alm_b}


# TODO: different gap penalties at the borders? -- strip border gaps
#       (related) benefit for longer non-gaps?
# TODO: allow exclusion, force inclusion? (exclusion probably must be done
# skipping over or setting an extremily high weight, and inclusion just
# joining different sub paths)
# TODO: gap opening and extension are traditionally defined as negative
# numbers, should we try to follow it here? The scorers are also following
# the tradition....
# TODO: should allow scoring gaps only in seq2, or in both?
# TODO: should have bidirectional
# TODO: should gap penalties be allowed to be functions?
# TODO: allow to search for *all* paths?
# TODO: do we really need to pass `seq_a` and `seq_b` again?
def align(
    graph, nodes, seq_a, seq_b, k, gap_open=1.0, gap_ext=0.0, gap="-", n_paths=None
):
    """
    Return the `k` best alignments in terms of costs.
    The function with take a path as source and target nodes in a graph,
    score the best alignments in terms of transition weights and gap
    penalties, and return the `k` best alignments. As an NP hard problem,
    it is not mathmatically guaranteed that the best path will be found
    unless the search pass, defined by `n_paths`, includes all possible
    paths.

    Parameters
    ==========
    graph : networkx object
        A directed weighted graph for the alignment, as returned by
        compute_graph().
    nodes : list
        A list of two graph labels in the format `(source, target)`,
        indicating the alignment path source and target.
    seq_a : list
        A list of elements used as the first sequence (on the vertical
        border).
    seq_b : list
        A list of elements used as the second sequences (on the horizontal
        border).
    k : int
        The maximum number of best alignments to be returned.
    gap_open : number
        The penalty score for gap openings. The higher the value, the
        more penalty will be applied to each gap. Defaults to 1.0.
    gap_ext : number
        The penalty score for gap extension. The higher the value, the
        more penalty will be applied to each gap. Defaults to 0.0.
    gap : str
        A string used as gap symbol. Defaults to `"-"`.
    n_paths : int
        The number of alignment paths to be collected for scoring. If
        not provided, a default value will be calculated based on the
        length of the sequences and the complexity of the graph.

    Returns
    =======
    alms : list
        A sorted list of best alignments in the format `[alm, cost]`,
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
    # TODO: change to consider `k` (or just scale with K)
    if not n_paths:
        n_paths = min(len(seq_a), len(seq_b))

    paths = list(
        islice(
            nx.shortest_simple_paths(graph, nodes[0], nodes[1], weight="weight"),
            n_paths,
        )
    )

    # Iterate over the collected paths, collecting the corresponding
    # alignments and weights so we can sort after the loop
    alignments = []
    for path in paths:
        # Sum edge weights
        weight = sum(
            [graph.edges[edge[1], edge[0]]["weight"] for edge in _pairwise_iter(path)]
        )

        # Build sequential representation of the alignment alignment
        alignment = build_align(path, seq_a, seq_b, gap=gap)

        # Add weights from gap opening and extension in both sequences
        # NOTE: as `groupby` returns an iterator, two subsequent list
        # comprehensions are needed in order not to consume it immediately,
        # thus casting to a list.
        gap_seqs_a = [list(g) for k, g in groupby(alignment["a"])]
        gap_seqs_b = [list(g) for k, g in groupby(alignment["b"])]
        gap_seqs_a = [g for g in gap_seqs_a if g[0] in gap]
        gap_seqs_b = [g for g in gap_seqs_b if g[0] in gap]

        # Update scores
        # TODO: decide how to normalize, as this is a weight
        alignment["score_a"] = weight + (gap_open * len(gap_seqs_a))
        alignment["score_a"] += sum([gap_ext * len(gap_seq) for gap_seq in gap_seqs_a])

        alignment["score_b"] = weight + (gap_open * len(gap_seqs_b))
        alignment["score_b"] += sum([gap_ext * len(gap_seq) for gap_seq in gap_seqs_b])

        alignment["score"] = alignment["score_a"] + alignment["score_b"]

        alignments.append(alignment)

    # sort by weight
    # TODO: improve sorting and make more reproducible (for same score)
    alignments = sorted(alignments, key=lambda a: a["score"])

    return alignments[:k]


# TODO: treat kwargs, etc
# TODO: decide on n_paths -- pass more than `k`, but how much?
def kbest_align(seq_a, seq_b, k=1, gap="-", matrix=None, **kwargs):

    if not matrix:
        matrix = utils.identity_matrix([seq_a, seq_b])

    graph = compute_graph(seq_a, seq_b, matrix)

    dest = "%i:%i" % (len(seq_a), len(seq_b))
    alms = align(graph, ("0:0", dest), seq_a, seq_b, k, n_paths=k * 2)

    return alms
