"""
Main module with code for alignment methods.
"""

# Import Python standard libraries
from collections import defaultdict
import itertools

# Import 3rd-party libraries
import numpy as np

# Import other modules
import malign.nw as nw
import malign.kbest as kbest  # TODO: rename to yenksp
import malign.utils as utils

# TODO: move to its own file, for simmetry
# TODO: expand dumb_malign by adding random gaps, call this pad_align?
def dumb_malign(seqs, gap="-", **kwargs):
    # Obtain the longest sequence length
    max_length = max([len(seq) for seq in seqs])

    # Pad all sequences in `alm`
    alm = {"seqs": []}
    scores = []
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
        scores.append(1.0 - (num_pad / max_length))

    # Add overall score
    alm["score"] = np.mean(scores)

    # The `dumb` method will always return a single aligment, but we
    # still return a list for compatibility with other methods
    return [alm]


# TODO: is the `gap` even needed? only in case we have no scorer?
def _malign(seqs, matrix, pw_func, gap="-", **kwargs):

    # If no matrix is provided, build an identity one
    if not matrix:
        matrix = utils.identity_matrix(seqs, match=+1, gap=-1)

    k = kwargs.get("k", 1)

    domains = list(itertools.combinations(range(len(seqs)), 2))  # was `pairs`

    # get submatrices
    # TODO: needed when creating a Matrix class?
    # TODO: only compute if needed -- for pairwise it is ready
    sub_matrix = matrix.compute_submatrices(domains)

    # Run pairwise alignment on all pairs, collecting all potential
    # alignments for each sequence in `potential`, where they are already
    # grouped by length
    # TODO: implement adding and removing gaps for +1/-1 offsets? -- this is related
    #       with the "do as many as possible" (replacing `longest`) below
    potential = defaultdict(lambda: defaultdict(set))
    for idx_x, idx_y in domains:
        # Run pairwise alignment
        alms = pw_func(seqs[idx_x], seqs[idx_y], k=k, matrix=sub_matrix[idx_x, idx_y])

        # Add by length
        for alm in alms:
            potential[len(alm["a"])][idx_x].add(tuple(alm["a"]))
            potential[len(alm["b"])][idx_y].add(tuple(alm["b"]))

    # Before taking the product of all potential alignments with longest
    # lenght, we need to make sure that all sequences have such length,
    # as there might be cases where all alignments were shorter; in order to
    # do so, we get all potentials with the longest alignment and align
    # again, this time using the gaps
    # TODO: don't consider only the longest -- as long as it is possible, do it,
    #       and score normally; this should also help guarantee we get as many
    #       alignments as requested as long as it is possible
    longest = max(potential)
    has_longest = list(potential[longest])

    # TODO: what if the alignment sequence is longer than longest? decide
    # whether to include the new one, realigning the top (most likely, but
    # more complex to implmenet) or just drop
    # TODO: make into its own function, perhaps in `utils`
    for seq_idx in range(len(seqs)):
        if seq_idx not in has_longest:
            for long_idx in has_longest:
                # aling the current one with all long_idx aligned sequences
                for aligned in potential[longest][long_idx]:
                    # make sure we get the correct order
                    if seq_idx < long_idx:
                        alms = pw_func(
                            seqs[seq_idx],
                            list(aligned),
                            k=k,
                            matrix=sub_matrix[seq_idx, long_idx],
                        )
                        for alm in alms:
                            potential[longest][seq_idx].add(tuple(alm["a"]))
                    else:
                        alms = pw_func(
                            list(aligned),
                            seqs[seq_idx],
                            k=k,
                            matrix=sub_matrix[long_idx, seq_idx],
                        )
                        for alm in alms:
                            potential[longest][seq_idx].add(tuple(alm["b"]))

    # Build all potential alignments and score them
    cand = [potential[longest][idx] for idx in sorted(potential[longest])]
    alms = []
    for aligns in itertools.product(*cand):
        alm = {"seqs": []}

        for entry in aligns:
            alm["seqs"].append(entry)

        # compute alignment score using matrix and gap opening/ext
        alm["score"] = np.mean([matrix[corr] for corr in zip(*alm["seqs"])])

        alms.append(alm)

    # sort so that the best comes first
    alms = sorted(alms, reverse=True, key=lambda e: e["score"])

    return alms


def dumb_align(seq_a, seq_b, gap="-", **kwargs):
    """
    Perform pairwise alignment with the `dumb` method (for testing purposes).
    """

    # TODO: remove when possible
    seq_a, seq_b = list(seq_a), list(seq_b)

    # Pre-compute padding character from lengths
    num_pad = abs(len(seq_a) - len(seq_b))
    left_pad = int(num_pad / 2)
    right_pad = num_pad - left_pad
    left_pad = [gap] * left_pad
    right_pad = [gap] * right_pad

    # Compute alignment vectors and extend the right one with gaps (if needed)
    if len(seq_a) < len(seq_b):
        alm_a = [*left_pad, *seq_a, *right_pad]
        alm_b = seq_b[:]
    elif len(seq_b) < len(seq_a):
        alm_a = seq_a[:]
        alm_b = [*left_pad, *seq_b, *right_pad]
    else:
        alm_a = seq_a[:]
        alm_b = seq_b[:]

    # Compute alignment from number of gaps
    # TODO: run actual scoring -- method is dumb, not scoring
    score = 1.0 - (num_pad / max([len(seq_a), len(seq_b)]))

    return [
        {"a": alm_a, "b": alm_b, "score": score, "score_a": score, "score_b": score}
    ]


def nw_align(seq_a, seq_b, **kwargs):
    """
    Perform pairwise alignment with the `nw` method.
    """

    # Get arguments
    k = kwargs.get("k", 1)
    m = kwargs.get("matrix", None)
    gap = kwargs.get("gap", "-")

    # Perform alignment
    alms = nw.nw_align(seq_a, seq_b, gap=gap, k=k, matrix=m)

    return alms


# TODO: treat kwargs, etc
# TODO: decide on n_paths -- pass more than `k`, but how much?
def kbest_align(seq_a, seq_b, k=1, gap="-", matrix=None, **kwargs):

    if not matrix:
        matrix = utils.identity_matrix([seq_a, seq_b])  # , 2, -2)

    graph = kbest.compute_graph(seq_a, seq_b, matrix)

    dest = "%i:%i" % (len(seq_a), len(seq_b))
    alms = kbest.align(graph, ("0:0", dest), seq_a, seq_b, k, n_paths=k * 2)

    return alms


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
        seq_a = list(seq_a)
        seq_b = list(seq_b)

    # Get default parameters
    gap = kwargs.get("gap", "-")
    k = kwargs.get("k", 1)
    matrix = kwargs.get("matrix", None)
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
        alms = nw_align(seq_a, seq_b, k=k, gap=gap, matrix=matrix)
    elif method == "kbest":
        alms = kbest_align(seq_a, seq_b, k=k, gap=gap)
    else:
        alms = dumb_align(seq_a, seq_b, gap=gap)

    return alms


# TODO: have a partial function for the single best alignment?
def multi_align(seqs, method, **kwargs):
    # Make sure all sequences are lists, as it facilitates manipulation later.
    # Note that, while not recommended, this even allows to mix different iterable
    # types in `seqs`.
    seqs = [list(seq) for seq in seqs]

    # Get default parameters
    gap = kwargs.get("gap", "-")
    k = kwargs.get("k", 1)
    matrix = kwargs.get("matrix", None)

    # Validate parameters
    if not gap:
        raise ValueError("Gap symbol must be a non-empty string.")
    if k < 1:
        raise ValueError("At least one alignment must be returned.")
    if method not in ["dumb", "nw", "yenksp"]:
        raise ValueError("Invalid alignment method `%s`." % method)

    # Run alignment method
    if method == "nw":
        pairwise_func = nw_align
        alms = _malign(seqs, matrix, pw_func=pairwise_func, gap=gap, k=k)
        # alms = _malign(seqs, kwargs["matrix"], pw_func=nw_align, gap=gap)
    elif method == "yenksp":
        pairwise_func = kbest_align
        alms = _malign(seqs, matrix, pw_func=pairwise_func, gap=gap, k=k)
        # alms = _malign(seqs, kwargs["matrix"], pw_func=kbest_align, gap=gap)
    else:
        alms = dumb_malign(seqs, gap=gap)

    #    alms = _malign(seqs, kwargs["matrix"], pw_func=pairwise_func, gap=gap)

    return alms
