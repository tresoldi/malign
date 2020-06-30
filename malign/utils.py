"""
Utility data and functions for the library.
"""

# Import Python standard libraries
import itertools
from string import ascii_uppercase

import numpy as np

# TODO: Remove temporary DNA scorer holder in future versions
DNA_MATRIX = {
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


def _label_iter():
    for length in itertools.count(1):
        for chars in itertools.product(ascii_uppercase, repeat=length):
            yield "".join(chars)


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


# TODO: implement sorting
def print_malms(alms):
    for alm_idx, alm in enumerate(alms):
        for label, seq in zip(_label_iter(), alm["seqs"]):
            print(f"{alm_idx} {label} ({alm['score']:.2f}) : {seq}")


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


# TODO: multiple matrices
def combine_matrices(matrix_a, matrix_b):
    # Collect the full alphabet in matrices A and B; assumes a[0] is the
    # same domain as b[0]
    alphabet_a = set([c[0] for c in matrix_a] + [c[0] for c in matrix_b])
    alphabet_b = set([c[1] for c in matrix_a])
    alphabet_c = set([c[1] for c in matrix_b])

    # Make copies of the matrices, so we don't change them
    matrix_a = matrix_a.copy()
    matrix_b = matrix_b.copy()

    def _add_missing_values(alphabet1, alphabet2, matrix):
        # TODO: here using the mean value

        for p in itertools.product(alphabet1, alphabet2):
            # mean values for ? <-> b
            value_b = {}
            if p not in matrix:
                # compute the value if missing
                if p[1] not in value_b:
                    values = [
                        value for key, value in matrix_a.items() if key[1] == p[1]
                    ]
                    if not values:
                        value_b[p[1]] = -1.0
                    else:
                        value_b[p[1]] = np.mean(values)

                # Add value to matrix
                matrix[p] = value_b[p[1]]

    # Make sure alphabet_a is mapped to all values in Ma and Mb
    _add_missing_values(alphabet_a, alphabet_b, matrix_a)
    _add_missing_values(alphabet_a, alphabet_c, matrix_b)

    # Build the new matrix as cleverly as it can
    matrix = {}
    for key in itertools.product(alphabet_a, alphabet_b, alphabet_c):
        val_a = matrix_a[key[0], key[1]]
        val_b = matrix_b[key[0], key[2]]
        matrix[key] = (val_a + val_b) / 2.0

    return matrix
