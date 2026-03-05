"""UPGMA-guided progressive alignment with beam search for k-best results.

Uses a UPGMA guide tree built from pairwise alignment distances to determine
the merge order. At each merge step, profiles are combined via gap insertion
guided by representative pairwise alignments, with beam search to maintain
the top-k candidates.
"""

import itertools
from collections.abc import Callable, Hashable

import numpy as np
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform

from .alignment import Alignment
from .ndim_common import should_use_ndim
from .ndim_yenksp import ndim_yenksp_align
from .scoring_matrix import ScoringMatrix
from .utils import score_alignment, sort_alignments


def _pairwise_distances(
    seqs: list[list[Hashable]],
    matrix: ScoringMatrix,
    pw_func: Callable[..., list[Alignment]],
) -> tuple[np.ndarray, dict[tuple[int, int], list[Alignment]]]:
    """Compute pairwise alignment distances and cache alignments.

    Args:
        seqs: Sequences to align pairwise.
        matrix: Scoring matrix with sub-matrices for all pairs.
        pw_func: Pairwise alignment function (yenksp_align or nw_align).

    Returns:
        Tuple of (NxN distance matrix, cache of pairwise alignments).
    """
    n = len(seqs)
    domains = list(itertools.combinations(range(n), 2))
    sub_matrices = matrix.compute_submatrices(domains)

    scores = np.zeros((n, n))
    cache: dict[tuple[int, int], list[Alignment]] = {}

    for i, j in domains:
        alms = pw_func(seqs[i], seqs[j], k=1, matrix=sub_matrices[i, j])
        cache[i, j] = alms
        if alms and alms[0].score is not None:
            scores[i, j] = alms[0].score
            scores[j, i] = alms[0].score

    max_score = scores.max() if scores.max() > 0 else 1.0
    dist = max_score - scores
    np.fill_diagonal(dist, 0.0)

    return dist, cache


def _expand_profile(
    profile: list[list[Hashable]],
    aligned_rep: list[Hashable],
    gap: Hashable,
) -> list[list[Hashable]]:
    """Expand profile columns to match an aligned representative.

    Uses a two-pointer approach: walks the aligned representative alongside
    the profile's representative (first sequence), inserting gap columns
    where the alignment introduced new gaps.

    Args:
        profile: List of sequences (all same length) forming the profile.
        aligned_rep: The aligned version of the profile's representative.
        gap: Gap symbol.

    Returns:
        Expanded profile with columns matching aligned_rep's length.
    """
    rep = profile[0]
    result = [[] for _ in profile]

    prof_col = 0
    for symbol in aligned_rep:
        if symbol == gap:
            # Alignment inserted a gap — add gap column to all profile seqs
            for seq_out in result:
                seq_out.append(gap)
        else:
            # Non-gap: consume next non-gap column from profile
            # Skip any profile-internal gap columns first (copy them through)
            while prof_col < len(rep) and rep[prof_col] == gap:
                for seq_idx, seq_out in enumerate(result):
                    seq_out.append(profile[seq_idx][prof_col])
                prof_col += 1

            # Now consume the matching non-gap column
            if prof_col < len(rep):
                for seq_idx, seq_out in enumerate(result):
                    seq_out.append(profile[seq_idx][prof_col])
                prof_col += 1

    # Copy any trailing gap columns from the profile
    while prof_col < len(rep):
        for seq_idx, seq_out in enumerate(result):
            seq_out.append(profile[seq_idx][prof_col])
        prof_col += 1

    return result


def _merge_profiles(
    left_profile: list[list[Hashable]],
    right_profile: list[list[Hashable]],
    matrix: ScoringMatrix,
    pw_func: Callable[..., list[Alignment]],
    k: int,
    gap: Hashable,
) -> list[list[list[Hashable]]]:
    """Merge two profiles, returning up to k alternative merged profiles.

    Aligns the representatives (first sequences) of each profile using
    the pairwise function, then expands both profiles to match the
    aligned representatives.

    Args:
        left_profile: Left profile (list of aligned sequences).
        right_profile: Right profile (list of aligned sequences).
        matrix: Pairwise scoring matrix for the representatives.
        pw_func: Pairwise alignment function.
        k: Number of alternative alignments to generate.
        gap: Gap symbol.

    Returns:
        List of up to k merged profiles (each is a list of sequences).
    """
    left_rep = [s for s in left_profile[0] if s != gap]
    right_rep = [s for s in right_profile[0] if s != gap]

    alms = pw_func(left_rep, right_rep, k=k, matrix=matrix)

    merged = []
    for alm in alms:
        aligned_left_rep = list(alm.seqs[0])
        aligned_right_rep = list(alm.seqs[1])

        expanded_left = _expand_profile(left_profile, aligned_left_rep, gap)
        expanded_right = _expand_profile(right_profile, aligned_right_rep, gap)

        # Pad all sequences to the same length (profiles may differ due to
        # internal gaps being preserved during expansion)
        combined = expanded_left + expanded_right
        max_len = max(len(s) for s in combined)
        padded = [s + [gap] * (max_len - len(s)) for s in combined]

        merged.append(padded)

    return merged


def _score_profile(
    profile: list[list[Hashable]],
    indices: list[int],
    matrix: ScoringMatrix,
) -> float:
    """Sum of all pairwise alignment scores within a profile.

    Args:
        profile: List of aligned sequences forming the profile.
        indices: Original sequence indices corresponding to profile rows.
        matrix: Full scoring matrix.

    Returns:
        Total score across all pairwise comparisons.
    """
    total = 0.0
    pairs = list(itertools.combinations(range(len(profile)), 2))
    if not pairs:
        return 0.0

    pair_domains = list(itertools.combinations(indices, 2))
    sub_matrices = matrix.compute_submatrices(pair_domains)

    for (pi, pj), (idx_i, idx_j) in zip(pairs, pair_domains, strict=True):
        sub_mtx = sub_matrices[idx_i, idx_j]
        total += score_alignment([profile[pi], profile[pj]], sub_mtx)

    return total


def upgma_progressive_align(
    seqs: list[list[Hashable]],
    matrix: ScoringMatrix,
    pw_func: Callable[..., list[Alignment]],
    k: int = 1,
) -> list[Alignment]:
    """UPGMA-guided progressive alignment with beam search for k-best.

    Builds a guide tree using UPGMA (average linkage) from pairwise
    alignment distances, then merges profiles bottom-up. At each merge
    step, beam search retains the top-k candidates.

    Args:
        seqs: List of sequences to align.
        matrix: Scoring matrix for all sequences.
        pw_func: Pairwise alignment function (yenksp_align or nw_align).
        k: Number of best alignments to return (beam width).

    Returns:
        Sorted list of up to k Alignment objects.
    """
    n = len(seqs)
    gap = matrix.gap

    if n < 2:
        return [Alignment(seqs=[list(s) for s in seqs], score=0.0)]

    if n == 2:
        return pw_func(seqs[0], seqs[1], k=k, matrix=matrix)

    # Step 1: Pairwise distances
    dist_matrix, _cache = _pairwise_distances(seqs, matrix, pw_func)

    # Step 2: UPGMA guide tree
    condensed = squareform(dist_matrix, checks=False)
    link = linkage(condensed, method="average")

    # Step 3: Bottom-up traversal with beam search
    # Each node stores a beam: list of (profile, indices) tuples
    # profile = list of aligned sequences, indices = original seq indices
    beams: dict[int, list[tuple[list[list[Hashable]], list[int]]]] = {}

    # Initialize leaves
    for i in range(n):
        beams[i] = [([list(seqs[i])], [i])]

    # Process merges from linkage matrix
    for merge_idx in range(len(link)):
        left_id = int(link[merge_idx, 0])
        right_id = int(link[merge_idx, 1])
        new_id = n + merge_idx

        left_beam = beams[left_id]
        right_beam = beams[right_id]

        # Collect all merged indices for this subtree
        sample_left_indices = left_beam[0][1]
        sample_right_indices = right_beam[0][1]
        merged_indices = sample_left_indices + sample_right_indices
        subtree_n = len(merged_indices)

        # Check if we can use true N-dim alignment for small subtrees
        if subtree_n <= 4 and should_use_ndim(
            subtree_n,
            [len(seqs[i]) for i in merged_indices],
            "yenksp",
        ):
            # Use true N-dim alignment on raw sequences
            subtree_seqs = [seqs[i] for i in merged_indices]
            sub_domains = list(itertools.combinations(merged_indices, 2))
            sub_matrices = matrix.compute_submatrices(sub_domains)
            # Build a local matrix for the subtree
            local_scores: dict[tuple[Hashable, ...], float] = {}
            for key, val in matrix.scores.items():
                if key is not None:
                    local_scores[key] = val

            try:
                ndim_alms = ndim_yenksp_align(subtree_seqs, k=k, matrix=matrix)
                new_beam = []
                for alm in ndim_alms:
                    profile = [list(s) for s in alm.seqs]
                    new_beam.append((profile, merged_indices))
                if new_beam:
                    beams[new_id] = new_beam[:k]
                    # Clean up
                    del beams[left_id]
                    del beams[right_id]
                    continue
            except Exception:
                pass  # Fall through to profile merging

        # Profile merging with beam search
        candidates: list[tuple[float, list[list[Hashable]], list[int]]] = []

        for left_profile, left_indices in left_beam:
            for right_profile, right_indices in right_beam:
                # Get pairwise sub-matrix for representatives
                li = left_indices[0]
                ri = right_indices[0]
                pair = (min(li, ri), max(li, ri))
                domains_list = [pair]
                sub_matrices = matrix.compute_submatrices(domains_list)
                sub_mtx = sub_matrices[pair]

                all_indices = left_indices + right_indices

                merged_profiles = _merge_profiles(
                    left_profile,
                    right_profile,
                    sub_mtx,
                    pw_func,
                    k=k,
                    gap=gap,
                )

                for profile in merged_profiles:
                    score = _score_profile(profile, all_indices, matrix)
                    candidates.append((score, profile, all_indices))

        # Keep top k by score
        candidates.sort(key=lambda x: x[0], reverse=True)
        beams[new_id] = [(prof, idx) for _, prof, idx in candidates[:k]]

        # Clean up
        del beams[left_id]
        del beams[right_id]

    # Step 4: Score final alignments with full N-domain matrix
    root_id = n + len(link) - 1
    root_beam = beams[root_id]

    results = []
    for profile, indices in root_beam:
        # Reorder profile to match original sequence order
        order = sorted(range(len(indices)), key=lambda i: indices[i])
        ordered_profile = [profile[i] for i in order]

        # Remove all-gap columns
        aln_len = len(ordered_profile[0])
        keep_cols = []
        for col in range(aln_len):
            column = tuple(ordered_profile[seq_idx][col] for seq_idx in range(len(ordered_profile)))
            if not all(s == gap for s in column):
                keep_cols.append(col)

        if keep_cols:
            cleaned = [[seq[col] for col in keep_cols] for seq in ordered_profile]
        else:
            cleaned = ordered_profile

        score = score_alignment(cleaned, matrix)
        results.append(Alignment(seqs=cleaned, score=score))

    return sort_alignments(results)
