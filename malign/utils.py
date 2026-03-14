"""Utility data and functions for the library."""

import itertools
from collections.abc import Hashable, Sequence
from string import ascii_uppercase

from tabulate import tabulate

from malign.scoring_matrix import ScoringMatrix

from .alignment import Alignment


def pad_sequence(seq, gap):
    """Prepend a gap symbol to a sequence."""
    return [gap, *list(seq)]


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


def score_alignment(
    seqs: Sequence[Sequence[Hashable]],
    scorer,
    *,
    gap="-",
    gap_ext=0.0,
    gap_open=0.0,
    normalize=False,
) -> float:
    """Compute the score of an alignment according to a scoring matrix.

    Args:
        seqs: Aligned sequences (all same length).
        scorer: Scoring matrix for looking up symbol-tuple scores.
        gap: Gap symbol. Defaults to "-".
        gap_ext: Gap extension penalty. Defaults to 0.0.
        gap_open: Gap opening penalty. Defaults to 0.0.
        normalize: Whether to normalize by alignment length. Defaults to False.

    Returns:
        The computed alignment score.
    """

    site_score = sum(scorer[corr] for corr in zip(*seqs, strict=False))

    grouped: list[list[list[Hashable]]] = [
        [list(g) for _, g in itertools.groupby(seq)] for seq in seqs
    ]
    gap_runs: list[list[int]] = [[len(g) for g in groups if g[0] == gap] for groups in grouped]
    gap_runs = [runs for runs in gap_runs if runs]

    seq_penalty = sum(sum(runs) * gap_ext for runs in gap_runs) + (len(gap_runs) * gap_open)
    total_score = site_score + seq_penalty

    if normalize and seqs and seqs[0]:
        return total_score / len(seqs[0])

    return total_score


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
