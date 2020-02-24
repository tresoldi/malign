from itertools import product

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
