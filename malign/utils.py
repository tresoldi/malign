"""
Utility data and functions for the library.
"""

# Import Python standard libraries
import itertools
from string import ascii_uppercase

# Import 3rd party tools
from tabulate import tabulate

# Import from the package
import malign.matrix as matrix

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


def sort_alignments(alms):
    """
    Sorts a list of alignments in-place.

    While the logic is simple, the usage of a single function guarantees that no
    competing implementations are design for each alignment method.
    """

    return sorted(alms, reverse=True, key=lambda e: (e["score"], e["seqs"]))


# TODO: gap extension as a function?
# TODO: in this case, we don't expect full gap vectors (that are really only
#       used for scoring), which here should be heavily penalized (or make
#       sure they are never colleted at all)
def score_alignment(seqs, scorer, **kwargs):
    """
    Returns the score of an alignment according to a matrix.
    """

    # Get parameters
    gap = kwargs.get("gap", "-")
    gap_ext = kwargs.get("gap_ext", -1)
    gap_open = kwargs.get("gap_open", -1.5)

    # Collect the scores for pure alignment sites
    site_score = [scorer[corr] for corr in zip(*seqs)]

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


# TODO: allow customizations
def tabulate_alms(alms):
    """
    Return a tabulated textual representation of alignments.
    """

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


# TODO: do sub-matrices and matrices at the same pass?
# TODO: add mismatch
def identity_matrix(seqs, **kwargs):
    """
    Build an identity matrix from a list of sequences.

    The function assumes, as expected, that the domains of all sequences are
    equal, even if symbols don't appear in all sequences.
    """

    # Collect scores
    match = kwargs.get("match", 1)
    gap = kwargs.get("gap", -1)
    mismatch = kwargs.get("mismatch", gap)
    gap_symbol = kwargs.get("gap_symbol", "-")

    # Collect alphabet and build key space
    alphabet = list(set(list(itertools.chain.from_iterable(seqs)) + [gap_symbol]))

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

    # Build matrix and return
    return matrix.ScoringMatrix(scores)
