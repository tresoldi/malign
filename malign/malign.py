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
import malign.dumb as dumb
import malign.yenksp as yenksp
import malign.utils as utils


# TODO: is the `gap` even needed? only in case we have no scorer?
def _malign(seqs, matrix, pw_func, gap="-", **kwargs):

    k = kwargs.get("k", 1)

    # Build a list of all paired domains
    domains = list(itertools.combinations(range(len(seqs)), 2))

    # Compute all the submatrices; while the same scores could be accessed directly
    # via the `matrix` (by using `None`s in the positions where the domain does not
    # apply) it makes debugging easier to pretend these are pairwise comparisons
    # with their own, non-multiwise aware, scorers
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

        # Store in `potential` by length
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

    # TODO: cut following `k`?

    return alms


# TODO: have a partial function for the single best alignment?
# TODO: gap opening/gap extension
def multi_align(seqs, method, **kwargs):
    # Make sure all sequences are lists, as it facilitates manipulation later.
    # Note that, while not recommended, this even allows to mix different iterable
    # types in `seqs`.
    seqs = [list(seq) for seq in seqs]

    # Get default parameters
    gap = kwargs.get("gap", "-")
    k = kwargs.get("k", 1)
    matrix = kwargs.get("matrix", None)

    # If no matrix was provided, build an identity one
    if not matrix:
        matrix = utils.identity_matrix(seqs, match=+1, gap=-1)

    # Validate parameters
    if not gap:
        raise ValueError("Gap symbol must be a non-empty string.")
    if k < 1:
        raise ValueError("At least one alignment must be returned.")
    if method not in ["dumb", "nw", "yenksp"]:
        raise ValueError("Invalid alignment method `%s`." % method)

    # Run alignment method; note that the `dumb` method does not rely in expansion
    # from pairwise alingments with `_malign` as others
    if method == "dumb":
        alms = dumb.dumb_malign(seqs, gap=gap)
    else:
        if method == "nw":
            pairwise_func = nw.nw_align
        elif method == "yenksp":
            # NOTE: This function also takes care of building graph, etc.
            pairwise_func = yenksp.yenksp_align

        alms = _malign(seqs, matrix, pw_func=pairwise_func, gap=gap, k=k)

    return alms
