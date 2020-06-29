# Import other modules
import malign.nw as nw
import malign.kbest as kbest

# TODO: score from 0 to 1 or from 1 to 0?


def dumb_malign(seqs, gap="-", **kwargs):
    # Obtain the longest sequence length
    max_length = max([len(seq) for seq in seqs])

    # Pad all sequences in `alm`
    alm = []
    for seq in seqs:
        # Computer lengths and bads
        num_pad = max_length - len(seq)
        left_pad_len = int(num_pad / 2)
        right_pad_len = num_pad - left_pad_len
        left_pad = [gap] * left_pad_len
        right_pad = [gap] * right_pad_len

        # Append the padded sequence and the score, here computed from the
        # number of gaps
        alm.append(
            {
                "seq": [*left_pad, *list(seq), *right_pad],
                "score": 1.0 - (num_pad / max_length),
            }
        )

    # The `dumb` method will always return a single aligment, but we
    # still return a list for compatibility with other methods
    return [alm]


# TODO: place gaps equally to borders
def dumb_align(seq_a, seq_b, gap="-", **kwargs):
    """
    Perform pairwise alignment with the `dumb` method (for testing purposes).
    """

    # Pre-compute padding character from lengths
    num_pad = abs(len(seq_a) - len(seq_b))
    left_pad = int(num_pad / 2)
    right_pad = num_pad - left_pad
    left_pad = [gap] * left_pad
    right_pad = [gap] * right_pad

    # Compute alignment vectors and extend the right one with gaps (if needed)
    # TODO: what if sequence is already a list by definition, from align()?
    alm_a = [token for token in seq_a]
    alm_b = [token for token in seq_b]

    if len(seq_a) < len(seq_b):
        alm_a = [*left_pad, *alm_a, *right_pad]
    else:
        alm_b = [*left_pad, *alm_b, *right_pad]

    # Compute alignment from number of gaps
    score = 1.0 - (num_pad / max([len(seq_a), len(seq_b)]))

    return [
        {"a": alm_a, "b": alm_b, "score": score, "score_a": score, "score_b": score}
    ]


# TODO: pass scorer
def nw_align(seq_a, seq_b, gap="-", **kwargs):
    """
    Perform pairwise alignment with the `nw` method.
    """

    alms = nw.nw_align(seq_a, seq_b, gap=gap, k=kwargs.get("k", 1))

    return alms


# TODO: treat kwargs, etc
# TODO: decide on n_paths
def kbest_align(seq_a, seq_b, k=1, gap="-", scorer=None, **kwargs):

    if not scorer:
        scorer = kbest.fill_scorer(set(seq_a), set(seq_b))

    graph = kbest.compute_graph(seq_a, seq_b, scorer)

    dest = "%i:%i" % (len(seq_a), len(seq_b))
    alms = kbest.align(graph, ("0:0", dest), seq_a, seq_b, k, n_paths=k * 2)

    return alms


# TODO: normalize/scale score
# TODO: have a partial function for the single best alignment?
def pw_align(seq_a, seq_b, **kwargs):
    """
    Return a sorted list of pairwise alignments with scores.

    The function takes two sequences `seq_a` and `seq_b` and returns
    a sorted list of the top `k` best alignments. An alignment is
    a dictionary with fields:

        - `a`, with the alignment of the first sequence
        - `b`, with the alignment of the second sequence
        - `score`, with the alignment overall score
        - `score_a`, with the alignment scored in terms of A to B
        - `score_b`, with the alignment scored in terms of B to A

    Parameters
    ==========
    seq_a : list
        The first sequence to be aligned.
    seq_b : list
        The second sequence to be aligned.
    k : int
        The maximum number of best alignments to be returned. Note that the
        returned value will be a list even if `k` is equal to one (that is,
        the single best alignment). Defaults to one.
    method : str
        The alignment method to be used; choices are are `"dumb"` (position
        based, intended for development and prototyping and always
        returning a single alignment), `"nw"` (asymmetric
        Needleman–Wunsch), and `kbest` (asymmetric graph k-best path).
    gap : str
        Gap symbol. Defaults to `"-"`.

    Returns
    =======
    alms : list
        A list of dictionaries with the top `k` alignments, as described
        above.
    """

    # TODO: remove need for this in future
    if isinstance(seq_a, str):
        seq_a = [c for c in seq_a]
        seq_b = [c for c in seq_b]

    # Get default parameters
    gap = kwargs.get("gap", "-")
    k = kwargs.get("k", 1)
    method = kwargs.get("method", "nw")

    # Validate parameters
    if not gap:
        raise ValueError("Gap symbol must be a non-empty string.")
    if k < 1:
        raise ValueError("At least one alignment must be returned.")
    if method not in ["dumb", "nw", "kbest"]:
        raise ValueError("Invalid alignment method `%s`." % method)

    # Run alignment method
    if method == "nw":
        alms = nw_align(seq_a, seq_b, k=k, gap=gap)
    elif method == "kbest":
        alms = kbest_align(seq_a, seq_b, k=k, gap=gap)
    else:
        alms = dumb_align(seq_a, seq_b, gap=gap)

    return alms


# TODO: add all arguments
def multi_align(seqs, method, **kwargs):
    # Get default parameters
    gap = kwargs.get("gap", "-")
    k = kwargs.get("k", 1)

    # Validate parameters
    if not gap:
        raise ValueError("Gap symbol must be a non-empty string.")
    if k < 1:
        raise ValueError("At least one alignment must be returned.")
    if method not in ["dumb", "nw", "kbest"]:
        raise ValueError("Invalid alignment method `%s`." % method)

    # Run alignment method
    if method == "nw":
        alms = None
    elif method == "kbest":
        alms = None
    else:
        alms = dumb_malign(seqs, gap=gap)

    return alms
