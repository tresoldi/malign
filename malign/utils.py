"""Utility data and functions for the library."""

import itertools
from collections import Counter
from collections.abc import Hashable, Sequence
from string import ascii_uppercase

from tabulate import tabulate

from malign.scoring_matrix import ScoringMatrix

from .alignment import Alignment

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


def sort_alignments(alms: list[Alignment]) -> list[Alignment]:
    """Sort a list of alignments by score (descending).

    Args:
        alms: List of alignments to sort.

    Returns:
        Sorted list of alignments (highest score first).
    """
    return sorted(alms, reverse=True, key=lambda e: (e.score, tuple(e.seqs)))


def score_alignment(seqs: Sequence[Sequence[Hashable]], scorer, **kwargs) -> float:
    """Compute the score of an alignment according to a scoring matrix.

    Args:
        seqs: Aligned sequences (all same length).
        scorer: Scoring matrix for looking up symbol-tuple scores.
        **kwargs: Optional gap, gap_ext, gap_open parameters.

    Returns:
        The computed alignment score.
    """
    gap = kwargs.get("gap", "-")
    gap_ext = kwargs.get("gap_ext", -1)
    gap_open = kwargs.get("gap_open", -1)

    site_score = sum(scorer[corr] for corr in zip(*seqs, strict=False))

    gap_seqs = [[list(g) for k, g in itertools.groupby(seq)] for seq in seqs]
    gap_seqs = [[len(g) for g in gap_seq if g[0] == gap] for gap_seq in gap_seqs]
    gap_seqs = [gap_seq for gap_seq in gap_seqs if gap_seq]

    seq_penalty = sum(sum(gap_seq) * gap_ext for gap_seq in gap_seqs) + (len(gap_seqs) * gap_open)

    return (site_score + seq_penalty) / len(seqs[0])


def tabulate_alms(alms) -> str:
    """Return a tabulated textual representation of alignments.

    Args:
        alms: List of Alignment objects to display.

    Returns:
        Formatted table string.
    """

    def _label_iter():
        for length in itertools.count(1):
            for chars in itertools.product(ascii_uppercase, repeat=length):
                yield "".join(chars)

    alm_len = max(len(alm.seqs[0]) for alm in alms)
    headers = ["Idx", "Seq", "Score"] + [f"#{i}" for i in range(alm_len)]
    colalign = tuple(["left", "left", "decimal"] + ["center"] * alm_len)
    table = []
    for alm_idx, alm in enumerate(alms):
        for label, seq in zip(_label_iter(), alm.seqs, strict=False):
            table.append([alm_idx, label, f"{alm.score:.2f}", *list(seq)])

        table.append(["" for _ in range(3 + alm_len)])

    table = table[:-1]

    return tabulate(table, headers=headers, colalign=colalign, tablefmt="github")


def identity_matrix(seqs: Sequence[Sequence[Hashable]], **kwargs) -> ScoringMatrix:
    """Build an identity matrix from a list of sequences.

    Args:
        seqs: List of sequences whose alphabets define the domains.
        **kwargs: Optional match, mismatch, gap, gap_score parameters.

    Returns:
        A ScoringMatrix with identity-style scoring.
    """
    match_score = kwargs.get("match", 1.0)
    gap_score = kwargs.get("gap_score", -1.0)
    gap = kwargs.get("gap", "-")
    mismatch_score = kwargs.get("mismatch", gap_score * 0.9)

    alphabet = list({*list(itertools.chain.from_iterable(seqs)), gap})

    scores = {}
    for key in itertools.product(alphabet, repeat=len(seqs)):
        counter = Counter(key)
        most_common = counter.most_common(1)[0]

        if most_common[0] == gap:
            scores[key] = gap_score
        else:
            scores[key] = most_common[1] ** (1 + match_score)

    domains = list(itertools.combinations(range(len(seqs)), 2))
    for domain in domains:
        for symbols in itertools.product(alphabet, alphabet):
            symbol_iter = iter(symbols)
            key = tuple(next(symbol_iter) if idx in domain else None for idx in range(len(seqs)))

            if gap in symbols:
                scores[key] = gap_score
            elif symbols[0] != symbols[1]:
                scores[key] = mismatch_score
            else:
                scores[key] = match_score

    return ScoringMatrix(scores)
