from collections import defaultdict
import itertools

import numpy as np

# Import other modules
import malign.nw as nw
import malign.kbest as kbest

# Import other modules
import malign.utils as utils


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


# TODO: move to utils later?
# TODO: decide on mean/median/top, second highest value
def compute_submatrices(matrix, pairs):
    # Collect/compute submatrices
    matrix_values = defaultdict(lambda: defaultdict(list))
    for pair in pairs:
        x, y = pair
        for key, value in matrix.items():
            matrix_values[x, y][key[x], key[y]].append(value)

    sub_matrix = {}
    for pair, data in matrix_values.items():
        sub_matrix[pair] = {
            key: np.percentile(values, 80) for key, values in data.items()
        }

    return sub_matrix


def _malign(seqs, matrix, pw_func, gap="-", **kwargs):

    k = kwargs.get("k", 1)

    domains = list(itertools.combinations(range(len(seqs)), 2))  # was `pairs`

    # get submatrices
    # TODO: needed when creating a Matrix class?
    # TODO: it is easier to create a new submatrix which can be passed around as
    #       "the original" than to deal with the various `None`s
    sub_matrix = matrix.compute_submatrices(domains)

    # Run pairwise alignment on all pairs, collecting all potential
    # alignments for each sequence in `potential`, where they are already
    # grouped by length
    # TODO: implement adding and removing gaps for +1/-1 offsets?
    potential = defaultdict(lambda: defaultdict(set))
    for pair in domains:  # TODO: rename, as might not be a pair
        x, y = pair

        # Run pairwise alignment
        alms = pw_func(seqs[x], seqs[y], k=k, matrix=sub_matrix[pair])

        print(pair, alms, len(alms[0]["a"]), len(alms[0]["b"]))
        print(sub_matrix[pair].scores)

        # Add by length
        for alm in alms:
            potential[len(alm["a"])][x].add(tuple(alm["a"]))
            potential[len(alm["b"])][y].add(tuple(alm["b"]))

    # Before taking the product of all potential alignments with longest
    # lenght, we need to make sure that all sequences have such length,
    # as there might be cases where all alignments were shorter; in order to
    # do so, we get all potentials with the longest alignment and align
    # again, this time using the gaps
    longest = max(potential)
    has_longest = list(potential[longest])

    # TODO: what if the alignment sequence is longer than longest? decide
    # whether to include the new one, realigning the top (most likely, but
    # more complex to implmenet) or just drop
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
    for a in itertools.product(*cand):
        alm = {"seqs": []}

        for entry in a:
            alm["seqs"].append(entry)

        # compute alignment score using matrix and gap opening/ext
        alm["score"] = np.mean([matrix[corr] for corr in zip(*alm["seqs"])])

        alms.append(alm)

    # sort so that the best comes first
    alms = sorted(alms, reverse=True, key=lambda e: e["score"])

    return alms


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

    alms = nw.nw_align(
        seq_a, seq_b, gap=gap, k=kwargs.get("k", 1), matrix=kwargs.get("matrix", None)
    )

    return alms


# TODO: treat kwargs, etc
# TODO: decide on n_paths
def kbest_align(seq_a, seq_b, k=1, gap="-", matrix=None, **kwargs):

    if not matrix:
        matrix = utils.fill_matrix(set(seq_a), set(seq_b))

    graph = kbest.compute_graph(seq_a, seq_b, matrix)

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


# TODO: add all arguments
def multi_align(seqs, method, **kwargs):
    # TODO: remove need for this in future
    if isinstance(seqs[0], str):
        seqs = [[c for c in seq] for seq in seqs]

    # Get default parameters
    gap = kwargs.get("gap", "-")
    k = kwargs.get("k", 1)

    # Validate parameters
    if not gap:
        raise ValueError("Gap symbol must be a non-empty string.")
    if k < 1:
        raise ValueError("At least one alignment must be returned.")
    if method not in ["dumb", "nw", "yenksp"]:
        raise ValueError("Invalid alignment method `%s`." % method)

    # Run alignment method
    # TODO: the cut for `k` should be performed by `_malign`, but in order
    # to collect more stuff, yenksp should probably have a higher value,
    # so we increase the serach space -- perhaps pass k*2 or k**2
    if method == "nw":
        alms = _malign(seqs, kwargs["matrix"], pw_func=nw_align, gap=gap)
    elif method == "yenksp":
        alms = _malign(seqs, kwargs["matrix"], pw_func=kbest_align, gap=gap)
    else:
        alms = dumb_malign(seqs, gap=gap)

    return alms
