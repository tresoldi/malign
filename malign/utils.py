"""
Utility data and functions for the library.
"""

# Import Python standard libraries
import itertools
from string import ascii_uppercase
from collections import Counter

import numpy as np
from tabulate import tabulate


from . import matrix

# TODO: Remove temporary DNA scorer holder in future versions
DNA_MATRIX = matrix.ScoringMatrix(
    {
        ("A", "A"): 10,
        ("A", "G"): -1,
        ("A", "C"): -3,
        ("A", "T"): -4,
        ("G", "A"): -1,
        ("G", "G"): 7,
        ("G", "C"): -5,
        ("G", "T"): -3,
        ("C", "A"): -3,
        ("C", "G"): -5,
        ("C", "C"): 9,
        ("C", "T"): 0,
        ("T", "A"): -4,
        ("T", "G"): -3,
        ("T", "C"): 0,
        ("T", "T"): 8,
        ("A", "-"): -5,
        ("G", "-"): -5,
        ("C", "-"): -5,
        ("T", "-"): -5,
        ("-", "A"): -5,
        ("-", "G"): -5,
        ("-", "C"): -5,
        ("-", "T"): -5,
    }
)


def _label_iter():
    for length in itertools.count(1):
        for chars in itertools.product(ascii_uppercase, repeat=length):
            yield "".join(chars)


def pairwise_iter(iterable):
    """
    Internal function for sequential pairwise iteration.

    The function follows the recipe in Python's itertools documentation.
    [https://docs.python.org/3/library/itertools.html]
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    """

    item_a, item_b = itertools.tee(iterable)
    next(item_a, None)

    return zip(item_a, item_b)


# TODO: gap extension as a function?
# TODO: in this case, we don't expect full gap vectors (that are really only
#       used for scoring), which here should be heavily penalized (or make
#       sure they are never colleted at all)
def score_alignment(seqs, matrix, **kwargs):
    """
    Returns the score of an alignment according to a matrix.
    """

    # Get parameters
    gap = kwargs.get("gap", "-")
    gap_ext = kwargs.get("gap_ext", -1)
    gap_open = kwargs.get("gap_open", -1.5)

    # Collect the scores for pure alignment sites
    site_score = [matrix[corr] for corr in zip(*seqs)]

    # Collect the gap sub-sequences for each sequence
    # 1st pass ->  [[['A'], ['T', 'T'], ['-'], ['C'], ['G', 'G'], ['A'], ['-', '-'] ...
    # 2nd pass ->  [[1, 2], [2]...] (if no gaps, `[ [], [] ]`)
    gap_seqs = [[list(g) for k, g in itertools.groupby(seq)] for seq in seqs]
    gap_seqs = [[len(g) for g in gap_seq if g[0] == gap] for gap_seq in gap_seqs]

    # Compute the penalty per sequence based on `gap_seqs`, and correct `site_score`
    seq_penalty = [
        sum(gap_seq) * gap_ext + len(gap_seq) * gap_open for gap_seq in gap_seqs
    ]

    return sum(site_score) + sum(seq_penalty)


# TODO: gap symbol (and check if not in alphabet)
# TODO: rename to indel?
# TODO: assume alphabet_b equal to alphabet_a if not provided?
def build_basic_matrix(alphabet_a, alphabet_b, match=1, mismatch=-1, gap=-1):
    gap_symbol = "-"

    # Initialize the scorer (with the default gap/gap equal to zero)
    # and add match and mismatch scores
    scorer = {(gap_symbol, gap_symbol): 0}
    for a, b in itertools.product(alphabet_a, alphabet_b):
        if a == b:
            scorer[a, b] = match
        else:
            scorer[a, b] = mismatch

    # Add A to gap and B to gap scores
    for a in alphabet_a:
        scorer[a, gap] = gap
    for b in alphabet_b:
        scorer[gap, b] = gap

    return scorer


# TODO: string alignment given negative numbers in score
def print_alms(alms):
    for idx, alm in enumerate(alms):
        print(
            "A %i (%.2f / %.2f):" % (idx, alm["score_a"], alm["score"]),
            " ".join(alm["a"]),
        )
        print(
            "B %i (%.2f / %.2f):" % (idx, alm["score_b"], alm["score"]),
            " ".join(alm["b"]),
        )


# TODO: implement sorting?
# TODO: allow customizations
def tabulate_alms(alms):
    if not alms:
        print("<No alignment>")
        return

    alm_len = len(alms[0]["seqs"][0])
    headers = ["Idx", "Seq", "Score"] + [f"#{i}" for i in range(alm_len)]
    colalign = tuple(["left", "left", "decimal"] + ["center"] * alm_len)
    table = []
    for alm_idx, alm in enumerate(alms):
        # Add a row for each sequence in the alignment
        for label, seq in zip(_label_iter(), alm["seqs"]):
            table.append([alm_idx, label, "%.2f" % alm["score"]] + list(seq))

        # Add blank row
        table.append(["" for i in range(3 + alm_len)])

    # Remove last empty blank row
    table = table[:-1]

    return tabulate(table, headers=headers, colalign=colalign, tablefmt="github")


# TODO: deal with potentially different gap symbols
# TODO: do sub-matrices and matrices at the same pass?
# TODO: add mismatch
def identity_matrix(seqs, **kwargs):
    # build a simple identity matrix, like in the voldemort example

    # Collect scores
    # TODO: rename to `xxx_scores`
    match = kwargs.get("match", 1)
    gap = kwargs.get("gap", -1)
    mismatch = kwargs.get("mismatch", gap)
    gap_symbol = kwargs.get("gap_symbol", "-")

    # Collect alphabet and build key space
    alphabet = list(set(list(itertools.chain.from_iterable(seqs)) + [gap_symbol]))
    space = [alphabet] * len(seqs)

    # Build pairwise sub-matrices first, using the identity logic
    scores = {}
    domains = list(itertools.combinations(range(len(seqs)), 2))
    for domain in domains:
        for symbols in itertools.product(alphabet, alphabet):
            symbol_iter = iter(symbols)
            key = tuple(
                [
                    next(symbol_iter) if idx in domain else None
                    for idx in range(len(seqs))
                ]
            )

            if gap_symbol in symbols:
                scores[key] = gap
            elif symbols[0] != symbols[1]:
                scores[key] = mismatch
            else:
                scores[key] = match

    # Build keys from all alphabets and fill scorer
    for key in itertools.product(*space):
        score = 0
        for symbol, count in Counter(key).items():
            # NOTE: remember that `gap` is a negative number
            if symbol == gap_symbol:
                score += count * gap
            elif count == 1:
                score += mismatch
            else:
                score += (count - 1) * match

    #        symbols = len(set(key))

    # Compute the score as a power of the match score, detracting gaps (note that
    # the gap score is supposed to be passed as a negative)
    #        score = (match) ** (len(seqs) - symbols)
    #        if key.count(gap_symbol):
    #            score += (key.count(gap_symbol) + 1) * gap
    #        scores[key] = score

    # Build matrix and return
    m = matrix.ScoringMatrix(scores)

    return m


# TODO: allow different gap symbol, also in align_graph
# TODO: allow to use median instead of mean (or even mode?)
# TODO: move to utility or scorer module
# TODO: fix comment defaults->kwargs
def fill_matrix(alpha_a, alpha_b, scorer=None, **kwargs):
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
    matches, mismatches, gaps = [], [], []
    for key, value in scorer.items():
        if key[0] == key[1]:
            matches.append(value)
        elif "-" in key:
            gaps.append(value)
        else:
            mismatches.append(value)

    match = kwargs.get("match", None)
    if not match:
        if matches:
            match = sum(matches) / len(matches)
        else:
            match = 1.0

    mismatch = kwargs.get("mismatch", None)
    if not mismatch:
        if mismatches:
            mismatch = sum(mismatches) / len(mismatches)
        else:
            mismatch = -1.0

    gap = kwargs.get("gap", None)
    if not gap:
        if gaps:
            gap = sum(gaps) / len(gaps)
        else:
            gap = -2.0

    # Build a product of all symbols and set the score if not found
    # TODO: check for gap symbol in alphabets and raise an error
    for key in itertools.product(alpha_a, alpha_b):
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


# TODO: multiple matrices and domains, now defaulting to (0,1) and (0,2)
# TODO: check if the gap symbols are the same
# TODO: allow to combine giving filenames? perhaps in matrix itself?
# TODO: can this be moved to inside ScoringMatrix? is it already there?
def combine_matrices(matrix_a, matrix_b):
    # Collect alphabets
    alphabets = [[], [], []]
    alphabets[0] = sorted(set(matrix_a.alphabets[0] + matrix_b.alphabets[0]))
    alphabets[1] = matrix_a.alphabets[1]
    alphabets[2] = matrix_b.alphabets[1]

    # we provide a single point ita/rus/grk
    sub_matrices = {(0, 1): matrix_a, (0, 2): matrix_b}

    # Initialize with an empty scorer and submatrices, besides alphabets
    # TODO: have a solution so `scores` is not necessary -- if not provided, computation
    #       of domain in `ScoringMatrix` will fail
    new_matrix = matrix.ScoringMatrix(
        scores={("-", "-", "-"): 0}, sub_matrices=sub_matrices, alphabets=alphabets
    )

    return new_matrix
