"""
Module with functions for the k-best pairwise alignment.
"""

# Import Python standard libraries
from itertools import product, islice, tee, groupby
import operator

# Import 3rd party libraries
import networkx as nx


# Define an internal auxiliary function.
def _pairwise(iterable):
    """
    Internal function for sequential pairwise iteration.
    The function follows the recipe in Python's itertools documentation.
    [https://docs.python.org/3/library/itertools.html]
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    """

    item_a, item_b = tee(iterable)
    next(item_a, None)

    return zip(item_a, item_b)


# TODO: allow different gap symbol, also in align_graph
# TODO: allow to use median instead of mean (or even mode?)
def fill_scorer(alpha_a, alpha_b, scorer=None, defaults=None):
    """
    Fill an incomplete scorer for multiple alignment.

    The function makes a copy of a dictionary `scorer`, if provided, and
    inserts missing key/value pairs. The values for missing
    match, mismatch, and gap alignments are either provided by means of
    the `defaults` argument or computed from existing values in
    `scorer`. In particular:

        - missing matching scores are set to the value of
          `default['match']` or, if unavailable, to the mean
          matching score in `scorer` (if no matching score is provided
          and it is impossible to compute the mean, the fallback score
          of 1.0 is used);
        - missing mismatching scores are set to the value of
          `default['mismatch']` or, if unavailable, to the mean
          mismatching score in `scorer` (if no mismatching score is provided
          and it is impossible to compute the mean, the fallback score
          of -1.0 is used);
        - missing gap scores are set to the value of
          `default['gap']` or, if unavailable, to the mean
          gap score in `scorer` (if no gap score is provided
          and it is impossible to compute the mean, the fallback score
          of -2.0 is used).

    The alignment of a gap with another gap, as in the top left corner of
    the Needleman-Wunsch matrix, is by default set or overridden to 0.0.

    Parameters
    ==========

    alpha_a : list
        A list with the symbols in alphabet A (the values in the first
        elements of the tuples serving as scorer keys).
    alpha_b : list
        A list with the symbols in alphabet B (the values in the second
        elements of the tuples serving as scorer keys).
    scorer : dict
        An optional alignment scorer as a dictionary with tuples of
        pairwise values as keys and numbers as scores. If not provided, a
        new scorer will be created.
    defaults : dict
        An optional dictionary with the default match, mismatch, and gap
        scores to be used. If not provided, scores will be inferred from
        `scorer`, if possible, or fall-back to the hard-coded values
        of 1.0, -1.0, and -2.0 respectively.

    Returns
    =======

    scorer : dict
        A copy of `scorer` or a new `scorer` with all possible alignments,
        including with gaps, specified.
    """

    # Either make a copy of the scorer, preserving the original,
    # or create a new one
    if not scorer:
        scorer = {}
    else:
        scorer = scorer.copy()

    # Set the score of ("-", "-") to the default
    scorer[("-", "-")] = 0.0

    # If `match`, `mismatch`, or `gap` were not provided, compute them from
    # their respective means, defaulting to 1.0, -1.0, and -2.0.
    # Note that, while this collects the match/mismatch/gap scores before
    # verifying if they are actually needed, considering that the first
    # fall-back has more logic involved (i.e., mean of values), it is still
    # better to perform this computation earlier, even if not needed,
    # than making the logic more complex to potentially save some cycles.
    # TODO: move to separate, util function
    defaults = defaults or {}
    matches, mismatches, gaps = [], [], []
    for key, value in scorer.items():
        if key[0] == key[1]:
            matches.append(value)
        elif "-" in key:
            gaps.append(value)
        else:
            mismatches.append(value)

    match = defaults.get("match", None)
    if not match:
        if matches:
            match = sum(matches) / len(matches)
        else:
            match = 1.0

    mismatch = defaults.get("mismatch", None)
    if not mismatch:
        if mismatches:
            mismatch = sum(mismatches) / len(mismatches)
        else:
            mismatch = -1.0

    gap = defaults.get("gap", None)
    if not gap:
        if gaps:
            gap = sum(gaps) / len(gaps)
        else:
            gap = -2.0

    # Build a product of all symbols and set the score if not found
    # TODO: check for gap symbol in alphabets and raise an error
    for key in product(alpha_a, alpha_b):
        if key not in scorer:
            if key[0] == key[1]:
                scorer[key] = match
            else:
                scorer[key] = mismatch

    # Set the gap score if not provided
    for symbol_a in alpha_a:
        if (symbol_a, "-") not in scorer:
            scorer[(symbol_a, "-")] = gap
    for symbol_b in alpha_b:
        if ("-", symbol_b) not in scorer:
            scorer[("-", symbol_b)] = gap

    return scorer


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
def compute_graph(seq_a, seq_b, scorer=None):
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

    As this alignment methods is designed as a building block for the
    alignment of sequences with different alphabets and with non symmetric
    transitions costs, the alignment and the scorer are intended to align
    the second sequence (`seq_b`) in terms of the first one (`seq_a`).
    Sequences can include gap symbols.

    Parameters
    ==========

    seq_a : list
        A list of hasheable elements to be aligned. The pairwise alignment
        is build in terms of this first sequence.
    seq_b : list
        A list of hasheable elements to be aligned.
    scorer: dict
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
    scorer = fill_scorer(set(seq_a), set(seq_b), scorer)
    max_score = max(scorer.values())

    # Add gaps to the beginning of both sequences, emulating the first row
    # and first column in NW.
    # NOTE: the "+" operation on lists here allows us to, in a single step,
    # add the necessary alignment gap and make an in-memory copy of both
    # sequences, preserving the original ones.
    seq_a = ["-"] + seq_a
    seq_b = ["-"] + seq_b

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
                ver_score = scorer[symbol_a, "-"]
            elif j == 0:
                dig_score = None
                hor_score = scorer["-", symbol_b]
                ver_score = None
            else:
                dig_score = scorer[symbol_a, symbol_b]
                hor_score = scorer["-", symbol_b]
                ver_score = scorer[symbol_a, "-"]

            # Add edges (and nodes automatically)
            if dig_score:
                graph.add_edge(
                    "%i:%i" % (i - 1, j - 1),
                    "%i:%i" % (i, j),
                    weight=max_score - dig_score,
                )
            if hor_score:
                graph.add_edge(
                    "%i:%i" % (i - 1, j),
                    "%i:%i" % (i, j),
                    weight=max_score - hor_score,
                )
            if ver_score:
                graph.add_edge(
                    "%i:%i" % (i, j - 1),
                    "%i:%i" % (i, j),
                    weight=max_score - ver_score,
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
    seq_a = seq_a[path[0][0] :]
    seq_b = seq_b[path[0][1] :]

    # Build the alignment sequences, adding (more) gaps if necessary
    alm_a = []
    alm_b = []
    for source, target in _pairwise(path):
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
def get_aligns(
    graph,
    nodes,
    seq_a,
    seq_b,
    k,
    gap_open=1.0,
    gap_ext=0.0,
    gap="-",
    n_paths=None,
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
    # to collect is difficult to determine beforehand, so wither the
    # user is allowed to provide it or we determine it with this simple
    # operation.
    if not n_paths:
        n_paths = min(len(seq_a), len(seq_b))

    paths = list(
        islice(
            nx.shortest_simple_paths(
                graph, nodes[0], nodes[1], weight="weight"
            ),
            n_paths,
        )
    )

    # Iterate over the collected paths, collecting the corresponding
    # alignments and weights so we can sort after the loop
    alignments = []
    for path in paths:
        # Sum edge weights
        weight = sum(
            [
                graph.edges[edge[1], edge[0]]["weight"]
                for edge in _pairwise(path)
            ]
        )

        # Build sequential representation of the alignment alignment
        alignment = build_align(path, seq_a, seq_b, gap=gap)

        # Add weights from gap opening and extension in both sequences
        # NOTE: as `groupby` returns an iterator, two subsequent list
        # comprehensions are needed in order not to consume it immediately,
        # thus casting to a list.
        for aligned_seq in alignment:
            # Get all the gap groups
            gap_seqs = [list(g) for k, g in groupby(aligned_seq)]
            gap_seqs = [g for g in gap_seqs if g[0] is gap]

            # Add weights due to gap opening and extension
            weight += gap_open * len(gap_seqs)
            weight += sum([gap_ext * len(gap_seq) for gap_seq in gap_seqs])

        alignments.append([alignment, weight])

    # sort by weight
    alignments = sorted(alignments, key=operator.itemgetter(1))

    return alignments[:k]
