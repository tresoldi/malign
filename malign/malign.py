# Import other modules
import malign.nw as nw

# TODO: score from 0 to 1 or from 1 to 0?

# TODO: make sure dumb_align returns a collection of a single alignment,
# as other methods do

# TODO: specify gap symbol
def dumb_align(seq_a, seq_b, args):
    """
    Perform pairwise alignment with the `dumb` method (for testing purposes).
    """

    # Pre-compute padding character from lengths
    num_pad = abs(len(seq_a) - len(seq_b))
    alm_pad = ["-"] * num_pad

    # Compute alignment vectors and extend the right one with gaps (if needed)
    # TODO: what if sequence is already a list by definition, from align()?
    alm_a = [token for token in seq_a]
    alm_b = [token for token in seq_b]

    if len(seq_a) < len(seq_b):
        alm_a += alm_pad
    else:
        alm_b += alm_pad

    # Compute alignment from number of gaps
    score = 1.0 - (num_pad / max([len(seq_a), len(seq_b)]))

    return alm_a, alm_b, score


# TODO: just returning the best for now, build later
def nw_align(seq_a, seq_b, args):
    """
    Perform pairwise alignment with the `nw` method.
    """

    seq_a = [c for c in seq_a]
    seq_b = [c for c in seq_b]
    alms = nw.nw_align(seq_a, seq_b)

    return alms[0]["a"], alms[0]["b"], alms[0]["score"]


# TODO: rename for pairwise and multiple
def align(seq_a, seq_b, method="nw", args=None):
    if method == "nw":
        alm_a, alm_b, score = nw_align(seq_a, seq_b, args)
    else:
        alm_a, alm_b, score = dumb_align(seq_a, seq_b, args)

    return alm_a, alm_b, score
