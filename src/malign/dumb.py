"""
Module for computing dumb (pure gap padding) alignments.
"""

# Import other modules
import malign.utils as utils
from .scoring_matrix import ScoringMatrix
from .alignment import Alignment

from typing import Sequence, Hashable, Optional


def dumb_malign(
    seqs: Sequence[Sequence[Hashable]],
    gap: Hashable = "-",
    matrix: Optional[ScoringMatrix] = None,
) -> Alignment:
    """
    Perform a *dumb* multiple alignment.

    This method is implemented for testing purposes, as it just pads gaps as necessary
    in order to return a single alignment.

    @param seqs:
    @param gap:
    @param matrix:
    @return:
    """

    # Get matrix, defaulting to an identity one
    if not matrix:
        matrix = utils.identity_matrix(seqs)

    # Obtain the longest sequence length
    max_length = max([len(seq) for seq in seqs])

    # Pad all sequences in `alm`
    ret_seqs = []
    for seq in seqs:
        # Computer lengths and bads
        num_pad = max_length - len(seq)
        left_pad_len = int(num_pad / 2)
        right_pad_len = num_pad - left_pad_len
        left_pad = [gap] * left_pad_len
        right_pad = [gap] * right_pad_len

        # Append the padded sequence and the score, here computed from the
        # number of gaps
        ret_seqs.append([*left_pad, *list(seq), *right_pad])

    # The `dumb` method will always return a single alignment, but we
    # still return a list for compatibility with other methods
    return Alignment(ret_seqs, utils.score_alignment(ret_seqs, matrix))
