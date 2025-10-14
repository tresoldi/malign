"""Module for computing dumb (pure gap padding) alignments."""

# Import Python standard libraries
from collections.abc import Hashable, Sequence

# Import other modules
from .alignment import Alignment
from .scoring_matrix import ScoringMatrix
from .utils import identity_matrix, score_alignment


# TODO: have a single signature accepting more than one sequence, and one for the alignment
# being the pairwise
def dumb_malign(
    seqs: Sequence[Sequence[Hashable]],
    matrix: ScoringMatrix | None = None,
) -> Alignment:
    """Perform a *dumb* multiple alignment.

    This method is implemented for testing purposes, as it just pads gaps as necessary
    in order to return a single alignment.

    @param seqs: List of sequences to be aligned.
    @param matrix: Scoring matrix.
    @return: The alignment for the two sequences.
    """

    # Get matrix, defaulting to an identity one
    if not matrix:
        matrix = identity_matrix(seqs)

    # Obtain the longest sequence length
    max_length = max([len(seq) for seq in seqs])

    # Pad all sequences in `alm`
    ret_seqs = []
    for seq in seqs:
        num_pad = max_length - len(seq)
        left_pad_len = int(num_pad / 2)
        right_pad_len = num_pad - left_pad_len
        left_pad = [matrix.gap] * left_pad_len
        right_pad = [matrix.gap] * right_pad_len

        # Append the padded sequence and the score, here computed from the
        # number of gaps
        ret_seqs.append([*left_pad, *list(seq), *right_pad])

    # The `dumb` method will always return a single alignment, but we
    # still return a list for compatibility with other methods
    return Alignment(ret_seqs, score_alignment(ret_seqs, matrix))
