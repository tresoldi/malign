"""Module for computing dumb (pure gap padding) alignments."""

from collections.abc import Hashable, Sequence

from .alignment import Alignment
from .scoring_matrix import ScoringMatrix
from .scoring_matrix import ScoringMatrix
from .utils import score_alignment


def dumb_malign(
    seqs: Sequence[Sequence[Hashable]],
    matrix: ScoringMatrix | None = None,
) -> Alignment:
    """Perform a dumb multiple alignment by padding with gaps.

    This baseline method just pads gaps as necessary to return a single
    alignment. Intended for testing purposes.

    Args:
        seqs: List of sequences to be aligned.
        matrix: Scoring matrix (default: identity matrix).

    Returns:
        A single Alignment with gap-padded sequences.
    """
    if not matrix:
        matrix = ScoringMatrix.from_sequences(seqs)

    max_length = max(len(seq) for seq in seqs)

    ret_seqs = []
    for seq in seqs:
        num_pad = max_length - len(seq)
        left_pad_len = int(num_pad / 2)
        right_pad_len = num_pad - left_pad_len
        left_pad = [matrix.gap] * left_pad_len
        right_pad = [matrix.gap] * right_pad_len

        ret_seqs.append([*left_pad, *list(seq), *right_pad])

    return Alignment(ret_seqs, score_alignment(ret_seqs, matrix))
