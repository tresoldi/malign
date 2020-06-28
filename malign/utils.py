"""
Utility data and functions for the library.
"""

# Import Python standard libraries
import itertools

# TODO: Remove temporary DNA scorer holder in future versions
DNA_SCORER = {
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

# TODO: gap symbol (and check if not in alphabet)
# TODO: rename to indel?
# TODO: assume alphabet_b equal to alphabet_a if not provided?
def build_basic_scorer(alphabet_a, alphabet_b, match=1, mismatch=-1, gap=-1):
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
            "A %i (%i / %i):" % (idx, alm["score_a"], alm["score"]), " ".join(alm["a"])
        )
        print(
            "B %i (%i / %i):" % (idx, alm["score_b"], alm["score"]), " ".join(alm["b"])
        )
