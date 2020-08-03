"""
Utility data and functions for the library.
"""

# Import Python standard libraries
from collections import Counter
from string import ascii_uppercase
import itertools

# Import 3rd party tools
from tabulate import tabulate

# Import from the package
from malign.scoring_matrix import ScoringMatrix

# TODO: Remove temporary DNA scorer holder in future versions
DNA_MATRIX = ScoringMatrix(
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
# TODO: add option to normalize score by length, which should probably be the default
#       given that we are going to compare and sort alignments of different lengths
# TODO: different gap penalties at the borders? -- strip border gaps
# TODO: (related) benefit for longer non-gaps?
def score_alignment(seqs, scorer, **kwargs):
    """
    Returns the score of an alignment according to a matrix.
    """

    # Get parameters
    gap = kwargs.get("gap", "-")
    gap_ext = kwargs.get("gap_ext", -1)
    gap_open = kwargs.get("gap_open", -1)

    # Collect the scores for pure alignment sites
    site_score = sum([scorer[corr] for corr in zip(*seqs)])

    # Collect the gap sub-sequences for each sequence
    # 1st pass ->  [[['A'], ['T', 'T'], ['-'], ['C'], ['G', 'G'], ['A'], ['-', '-'] ...
    # 2nd pass ->  [[1, 2], [2]...] (if no gaps, `[ [], [] ]`)
    # 3rd pass -> remove empty lists
    gap_seqs = [[list(g) for k, g in itertools.groupby(seq)] for seq in seqs]
    gap_seqs = [[len(g) for g in gap_seq if g[0] == gap] for gap_seq in gap_seqs]
    gap_seqs = [gap_seq for gap_seq in gap_seqs if gap_seq]

    # Compute the penalty per sequence based on `gap_seqs`, and correct `site_score`
    seq_penalty = sum([sum(gap_seq) * gap_ext for gap_seq in gap_seqs]) + (
        len(gap_seqs) * gap_open
    )

    # Correct by length, as we may be comparing alignments of different length
    return (site_score + seq_penalty) / len(seqs[0])


# TODO: allow customizations
def tabulate_alms(alms):
    """
    Return a tabulated textual representation of alignments.
    """

    # Internal function for generating representation labels
    def _label_iter():
        for length in itertools.count(1):
            for chars in itertools.product(ascii_uppercase, repeat=length):
                yield "".join(chars)

    alm_len = max([len(alm["seqs"][0]) for alm in alms])
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
def identity_matrix(seqs, **kwargs):
    """
    Build an identity matrix from a list of sequences.

    The function assumes, as expected, that the domains of all sequences are
    equal, even if symbols don't appear in all sequences.
    """

    # Collect scores
    match_score = kwargs.get("match", 1.0)
    gap_score = kwargs.get("gap_score", -1.0)
    gap = kwargs.get("gap", "-")
    #    mismatch_score = kwargs.get("mismatch", (match_score + gap_score) / 2.0)
    mismatch_score = kwargs.get("mismatch", gap_score * 0.9)  # TODO: just a bit below

    # Collect alphabet and build key space
    alphabet = list(set(list(itertools.chain.from_iterable(seqs)) + [gap]))

    # Build full scoring matrix
    scores = {}
    for key in itertools.product(alphabet, repeat=len(seqs)):
        counter = Counter(key)

        most_common = counter.most_common(1)[0]

        if most_common[0] == gap:
            scores[key] = gap_score
        else:
            # TODO: review, the power and +1 might be too much
            scores[key] = most_common[1] ** (1 + match_score)

    # Build pairwise sub-matrices first, using the identity logic
    # TODO: redo considering the full key from the number of seqs
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

            if gap in symbols:
                scores[key] = gap_score
            elif symbols[0] != symbols[1]:
                scores[key] = mismatch_score
            else:
                scores[key] = match_score

    # Build matrix and return
    return ScoringMatrix(scores)
