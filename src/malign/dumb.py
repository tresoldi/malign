"""Module for computing dumb (pure gap padding) alignments."""

# Import Python standard libraries
from typing import Hashable, Optional, Sequence

# Import other modules
from .alignment import Alignment
from .scoring_matrix import ScoringMatrix
from .utils import identity_matrix, score_alignment


def dumb_malign(
    seqs: Sequence[Sequence[Hashable]],
    matrix: Optional[ScoringMatrix] = None,
) -> Alignment:
    """
    Perform a *dumb* multiple alignment.

    This method is implemented for testing purposes and for getting the quickest
    alignment (where all sequences have then same length), as it just pads gaps as
    necessary in order to return a single alignment.

    @param seqs: List of sequences to be aligned.
    @param matrix: Scoring matrix.
    @return: The alignment for the two sequences.
    """

    # Get matrix, defaulting to an identity one; for the `dumb` alignment,
    # this is only used for scoring and for providing the gap symbol
    if not matrix:
        matrix = identity_matrix(seqs)

    # Obtain the longest sequence length
    max_length = max([len(seq) for seq in seqs])

    # Pad all sequences in `alm`
    ret_seqs = []
    for seq in seqs:
        # Compute left and right pad length for the sequence, then append
        # the padded sequence
        num_pad = max_length - len(seq)
        left_pad_len = int(num_pad / 2)

        left_pad = [matrix.gap] * left_pad_len
        right_pad = [matrix.gap] * (num_pad - left_pad_len)

        ret_seqs.append([*left_pad, *list(seq), *right_pad])

    # The `dumb` method will always return a single alignment, but we
    # still return a list for compatibility with other methods
    return Alignment(ret_seqs, score_alignment(ret_seqs, matrix))
