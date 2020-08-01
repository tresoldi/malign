"""
Module for computing dumb (pure gap padding) alignments.
"""

# Import other modules
import malign.utils as utils


def dumb_malign(seqs, gap="-", **kwargs):
    """
    Perform a *dumb* multiple alignment.

    This method is implemented for testing purposes, as it just pads gaps as necessary
    in order to return a single alignment.
    """

    # Get matrix, defaulting to an identity one
    matrix = kwargs.get("matrix", None)
    if not matrix:
        matrix = utils.identity_matrix(seqs)

    # Obtain the longest sequence length
    max_length = max([len(seq) for seq in seqs])

    # Pad all sequences in `alm`
    alm = {"seqs": []}
    for seq in seqs:
        # Computer lengths and bads
        num_pad = max_length - len(seq)
        left_pad_len = int(num_pad / 2)
        right_pad_len = num_pad - left_pad_len
        left_pad = [gap] * left_pad_len
        right_pad = [gap] * right_pad_len

        # Append the padded sequence and the score, here computed from the
        # number of gaps
        alm["seqs"].append([*left_pad, *list(seq), *right_pad])

    # Add overall score
    alm["score"] = utils.score_alignment(seqs, matrix)

    # The `dumb` method will always return a single aligment, but we
    # still return a list for compatibility with other methods
    return [alm]
