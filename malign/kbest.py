# Import Python standard libraries
from itertools import product

# Import 3rd party libraries
import networkx as nx

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
