"""
Main module with code for alignment methods.
"""

# Import Python standard libraries
from collections import defaultdict
import itertools

# Import other modules
import malign.anw as anw
import malign.dumb as dumb
import malign.yenksp as yenksp
import malign.utils as utils


def _build_candidates(potential_alms, matrix):
    """
    Internal function used by `_malign()`.

    This function takes a collection of individual sequence alignments, grouped by
    lengths in `potential_alms`, and a scoring `matrix`, and returns all potential
    alingments that can be constructed from it.

    Note that the function does not score the alignments.
    """

    # Build set for collecting alignments, num_seqs is computed only once
    alms = set()
    num_seqs = None

    # Build all potential alignments and score them
    cand = [potential_alms[idx] for idx in sorted(potential_alms)]
    for aligns in itertools.product(*cand):
        # Get the sequences and, if `num_seqs` has not been computed yet, cache it
        # (as it the same for all alignments)
        seqs = tuple(aligns)
        if not num_seqs:
            num_seqs = len(seqs)

        # The combination of alignments might results in full-gap vectors; here we
        # remove them and make sure to only add unique alignments
        full_gap = tuple(matrix.gap * num_seqs)
        vectors = list(zip(*seqs))
        if full_gap in vectors:
            vectors = [vector for vector in vectors if vector != full_gap]
            seqs = tuple(zip(*vectors))

        # Store the sequences, either in the original form or cleaned of full gaps
        alms.add(seqs)

    return alms


# pylint: disable=too-many-locals
def _collect_alignments(seqs, matrix, pw_func, **kwargs):
    """
    Internal function for multiwise alignment.

    This function is used to perform the multiwise alignment with different
    pairwise methods. It follows a procedure different from the more common
    UPGMA/NJ guiding trees, as it computes all potential pairwise alignments and
    expand it to full multiwise alignment for a single-scoring round.
    """

    # Collect additional arguments; `k` defaults to `None`, so that each pairwise
    # alignment function can decide its value
    gap = kwargs.get("gap", matrix.gap)
    k = kwargs.get("k", None)

    # Build a list of all paired domains
    # -> (0, 1), (0, 2), (1, 2)...
    domains = list(itertools.combinations(range(len(seqs)), 2))

    # Compute all the submatrices; while the same scores could be accessed directly
    # via the `matrix` (by using `None`s in the positions where the domain does not
    # apply) it makes debugging easier to "pretend" these are pairwise comparisons
    # with their own, non-multiwise aware, scorers.
    sub_matrix = matrix.compute_submatrices(domains)

    # Run pairwise alignment on all pairs, collecting all potential
    # alignments for each sequence in `potential`, where they are already
    # grouped by length. Note that the alignments must be converted to tuples in
    # order to be added to the `potential` sets (as lists are non-hasheable).
    potential = defaultdict(lambda: defaultdict(set))
    for idx_x, idx_y in domains:
        # Run pairwise alignment
        alms = pw_func(
            seqs[idx_x], seqs[idx_y], k=k, gap=gap, matrix=sub_matrix[idx_x, idx_y]
        )

        # Store in `potential` by length
        for alm in alms:
            potential[len(alm["seqs"][0])][idx_x].add(tuple(alm["seqs"][0]))
            potential[len(alm["seqs"][1])][idx_y].add(tuple(alm["seqs"][1]))

    # Starting from the minimum length (the maximum sequence length), fill all the
    # potential alignments
    min_length = max([len(seq) for seq in seqs])
    for length in sorted(potential):
        # Skip if less than minimum or already full
        if length <= min_length:
            continue
        if len(potential[length]) == len(seqs):
            continue

        # Get list of all indexes with the current length and list of those to compute;
        # `has_longest` and `to_compute` must be computed before the loop, as the
        # dictionary will change
        has_longest = list(potential[length])
        to_compute = [idx for idx in range(len(seqs)) if idx not in potential[length]]

        # Compute as necessary
        for seq_idx, long_idx in itertools.product(to_compute, has_longest):
            # aling the current one with all long_idx aligned sequences
            for aligned in potential[length][long_idx]:
                # Make sure we get the correct order; note that the calling is
                # more complex than it would be expected as we need to account
                # for the asymmetry in `sub_matrix` indexing
                if seq_idx < long_idx:
                    seq_a = seqs[seq_idx]
                    seq_b = list(aligned)
                    mtx = sub_matrix[seq_idx, long_idx]
                    alm_idx = 0  # seqs[seq_idx] is the first element
                else:
                    seq_a = list(aligned)
                    seq_b = seqs[seq_idx]
                    mtx = sub_matrix[long_idx, seq_idx]
                    alm_idx = 1  # seqs[seq_idx] is the second element

                # Align and add only those with the requested length
                for alm in pw_func(seq_a, seq_b, gap=gap, k=k, matrix=mtx):
                    if len(alm["seqs"][alm_idx]) == length:
                        potential[length][seq_idx].add(tuple(alm["seqs"][alm_idx]))

    # Build all candidate alignments, sort, and return
    alms = set()
    for length in potential:
        if len(potential[length]) == len(seqs):
            alms = alms.union(_build_candidates(potential[length], matrix))

    # Compute scores, sort and return
    alms = [
        {"seqs": seqs, "score": utils.score_alignment(seqs, matrix)} for seqs in alms
    ]

    return utils.sort_alignments(alms)


# TODO: gap opening/gap extension for scoring
def multi_align(seqs, method, **kwargs):
    """
    Compute multiple alignments for a list of sequences.

    The function will return a sorted list of the `k` best alignments, as long as they
    can be computed.

    Parameters
    ==========
    seqs : list of iterables
        A list of iterables (lists, tuples, or strings) representing the sequences to
        be aligned.
    method : str
        The method to be used for alignment computation. Currently supported methods are
        `"dumb"`, `"anw"` (for "asymmetric Needleman–Wunsch"), and `"yenksp"` (for
        the graph-alignment based on Yen's k-shortest paths algorithm).
    matrix :  dict or ScoringMatrix
        The matrix used for scoring the alignment. If provided, must match in length the
        number of sequences in `seq`. If not provided, an identity matrix will be created
        and used.
    k : int
        The maximum number of alignments to return. As there is no guarantee that the
        method being used or the sequences provided allow for `k` different alignments,
        the actual number might be less than `k`. Defaults to 1.
    gap : string
        The symbol to be used for gap representation. If `matrix` is a ScoringMatrix,
        it must match the gap symbol specified in it. Defaults to `"-"`.
    """

    # Make sure all sequences are lists, as it facilitates later manipulation.
    # Note that, while not recommended, this even allows to mix different iterable
    # types in `seqs`. The conversion also guarantees that a copy is made, for
    # immutability.
    seqs = [list(seq) for seq in seqs]

    # Get the user-provided matrix, or compute an identity one
    matrix = kwargs.get("matrix", None)
    if not matrix:
        matrix = utils.identity_matrix(seqs, match=+1, gap_score=-1)

    # Get parameters and validate them
    gap = kwargs.get("gap", "-")
    if not gap:
        raise ValueError("Gap symbol must be a non-empty string.")

    if not isinstance(matrix, dict):  # ScoringMatrix
        if gap != matrix.gap:
            raise ValueError("Different gap symbols.")

    k = kwargs.get("k", 1)
    if k < 1:
        raise ValueError("At least one alignment must be returned.")

    if method not in ["dumb", "anw", "yenksp"]:
        raise ValueError(f"Invalid alignment method `{method}`.")

    # Run alignment method; note that the `dumb` method does not rely in expansion
    # from pairwise alingments with `_malign` as others
    if method == "dumb":
        alms = dumb.dumb_malign(seqs, gap=gap)
    else:
        if method == "anw":
            pairwise_func = anw.nw_align
            pw_k = k
        elif method == "yenksp":
            # For `yenksp`, we will compute the twice the number of paths
            # requested, in order to get out of local minima
            pairwise_func = yenksp.yenksp_align
            pw_k = k ** 2

        alms = _collect_alignments(seqs, matrix, pw_func=pairwise_func, gap=gap, k=pw_k)

    return alms[:k]
