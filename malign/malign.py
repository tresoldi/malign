"""Main module with code for alignment methods."""

import itertools
from collections import defaultdict
from collections.abc import Callable, Hashable

from .alignment import Alignment
from .anw import nw_align
from .blocks import merge_alignment_blocks
from .dumb import dumb_malign
from .ndim_common import should_use_ndim
from .ndim_yenksp import ndim_yenksp_align
from .progressive import upgma_progressive_align
from .scoring_matrix import ScoringMatrix
from .utils import score_alignment, sort_alignments
from .yenksp import yenksp_align

_PW_FUNCS = {"anw": nw_align, "yenksp": yenksp_align}


def _build_candidates(
    potential_alms: dict[int, set[tuple[Hashable, ...]]],
    matrix: ScoringMatrix,
) -> set[tuple[tuple[Hashable, ...]]]:
    """Build all potential multiwise alignment candidates from pairwise results.

    Args:
        potential_alms: Dictionary with sequence indexes as keys and sets of
            collected potential alignments as values.
        matrix: The scoring matrix for all potential alignments.

    Returns:
        Set of all potential alignment combinations.
    """
    alms = set()
    num_seqs = None

    cand = [potential_alms[idx] for idx in sorted(potential_alms)]
    for aligns in itertools.product(*cand):
        seqs = tuple(aligns)
        if not num_seqs:
            num_seqs = len(seqs)

        full_gap = tuple([matrix.gap] * num_seqs)
        vectors = list(zip(*seqs, strict=False))
        if full_gap in vectors:
            vectors = [vector for vector in vectors if vector != full_gap]
            seqs = tuple(zip(*vectors, strict=False))

        alms.add(seqs)

    return alms


def _collect_alignments(
    seqs: list[list[Hashable]],
    matrix: ScoringMatrix,
    pw_func: Callable[..., list[Alignment]],
    k: int | None = None,
) -> list[Alignment]:
    """Perform multiwise alignment by expanding pairwise results.

    Args:
        seqs: List of sequences to be aligned.
        matrix: Scoring matrix (must match number of sequences).
        pw_func: Pairwise alignment function.
        k: Maximum number of alignments to return per pair.

    Returns:
        Sorted list of all collected Alignments.
    """
    domains = list(itertools.combinations(range(len(seqs)), 2))

    sub_matrix = matrix.compute_submatrices(domains)

    potential = defaultdict(lambda: defaultdict(set))
    for idx_x, idx_y in domains:
        alms = pw_func(seqs[idx_x], seqs[idx_y], k=k, matrix=sub_matrix[idx_x, idx_y])

        for alm in alms:
            seq_a, seq_b = tuple(alm.seqs[0]), tuple(alm.seqs[1])
            potential[len(seq_a)][idx_x].add(seq_a)
            potential[len(seq_b)][idx_y].add(seq_b)

    min_length = max(len(seq) for seq in seqs)
    for length in sorted(potential):
        if length <= min_length:
            continue
        if len(potential[length]) == len(seqs):
            continue

        has_longest = list(potential[length])
        to_compute = [idx for idx in range(len(seqs)) if idx not in potential[length]]

        for seq_idx, long_idx in itertools.product(to_compute, has_longest):
            for aligned in potential[length][long_idx]:
                if seq_idx < long_idx:
                    seq_a = seqs[seq_idx]
                    seq_b = list(aligned)
                    mtx = sub_matrix[seq_idx, long_idx]
                    alm_idx = 0
                else:
                    seq_a = list(aligned)
                    seq_b = seqs[seq_idx]
                    mtx = sub_matrix[long_idx, seq_idx]
                    alm_idx = 1

                for alm in pw_func(seq_a, seq_b, k=k, matrix=mtx):
                    if len(alm.seqs[alm_idx]) == length:
                        potential[length][seq_idx].add(tuple(alm.seqs[alm_idx]))

    alms = set()
    for length in potential:
        if len(potential[length]) == len(seqs):
            alms = alms.union(_build_candidates(potential[length], matrix))

    ret_alms = [Alignment(seqs, score_alignment(seqs, matrix)) for seqs in alms]

    return sort_alignments(ret_alms)


def align(
    sequences: list[Hashable],
    method: str = "anw",
    matrix: ScoringMatrix | None = None,
    k: int = 1,
    merge_blocks: bool = False,
    max_block_size: int = 2,
) -> list[Alignment]:
    """Compute multiple alignments for a list of sequences.

    Args:
        sequences: List of sequences to be aligned (strings or lists of symbols).
        method: Alignment method - "dumb", "anw", or "yenksp" (default: "anw").
        matrix: Scoring matrix. If not provided, an identity matrix is created.
        k: Maximum number of alignments to return (default: 1).
        merge_blocks: If True, merge complementary-gap block columns into
            compound symbols (tuples) before returning (default: False).
        max_block_size: Maximum number of columns in a detected block
            (default: 2). Only used when merge_blocks is True.

    Returns:
        Sorted list of the k-best multiple alignments.

    Raises:
        ValueError: If k < 1 or method is invalid.
    """
    seqs: list[list[Hashable]] = [list(seq) for seq in sequences]

    if not matrix:
        matrix = ScoringMatrix.from_sequences(seqs, match=+1, gap_score=-1)

    if k < 1:
        raise ValueError("At least one alignment must be returned.")

    if method not in ["dumb", "anw", "yenksp"]:
        raise ValueError(f"Invalid alignment method `{method}`.")

    if method == "dumb":
        alms = [dumb_malign(seqs, matrix=matrix)]
    else:
        pw_func = _PW_FUNCS[method]
        if len(seqs) == 2:
            alms = pw_func(seqs[0], seqs[1], k=k, matrix=matrix)
        elif should_use_ndim(len(seqs), [len(s) for s in seqs], method):
            alms = ndim_yenksp_align(seqs, k=k, matrix=matrix)
        else:
            alms = upgma_progressive_align(seqs, matrix, pw_func=pw_func, k=k)

    if merge_blocks:
        gap = matrix.gap if matrix else "-"
        alms = [merge_alignment_blocks(a, gap=gap, max_block_size=max_block_size) for a in alms]

    return alms
