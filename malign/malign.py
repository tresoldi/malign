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
