"""Matrix learning from cognate sets.

This module provides methods to learn scoring matrices from collections
of cognate sets (sequences believed to be related). The learned matrices
optimize alignment quality across all provided cognates.

Core Assumption: If cognate sets are correct, their alignment scores should
be collectively maximized under the optimal scoring matrix.
"""

from collections import Counter
from collections.abc import Hashable, Sequence

import numpy as np
from scipy.optimize import minimize

from .alignment import Alignment
from .malign import multi_align
from .scoring_matrix import ScoringMatrix


def learn_matrix(
    cognate_sets: list[list[Sequence[Hashable]]],
    method: str = "em",
    max_iter: int = 10,
    initial_matrix: ScoringMatrix | None = None,
    gap: Hashable = "-",
    **kwargs,
) -> ScoringMatrix:
    """Learn a scoring matrix from cognate sets.

    Given collections of related sequences (cognates), learns a scoring matrix
    that maximizes alignment quality across all sets.

    Args:
        cognate_sets: List of cognate sets, where each set is a list of sequences
            believed to be related (e.g., [["word1a", "word1b"], ["word2a", "word2b"]]).
        method: Learning method - "em" or "gradient_descent" (default: "em").
        max_iter: Maximum iterations for learning (default: 10).
        initial_matrix: Starting matrix (if None, creates identity matrix).
        gap: Gap symbol (default: "-").
        **kwargs: Additional method-specific parameters.

    Returns:
        Learned ScoringMatrix optimized for the provided cognates.

    Example:
        >>> cognate_sets = [
        ...     [["ACGT"], ["AGCT"]],
        ...     [["TGCA"], ["TGGA"]],
        ... ]
        >>> matrix = learn_matrix(cognate_sets, method="em", max_iter=5)

    Note:
        Phase 2 implements basic iteration without convergence criteria.
        Convergence checking will be added in Phase 3.
    """
    if method == "em":
        return _em_learning(cognate_sets, max_iter, initial_matrix, gap, **kwargs)
    elif method == "gradient_descent":
        return _gradient_descent_learning(cognate_sets, max_iter, initial_matrix, gap, **kwargs)
    else:
        raise ValueError(f"Unknown learning method: {method}. Use 'em' or 'gradient_descent'.")


def _initialize_matrix(
    cognate_sets: list[list[Sequence[Hashable]]], gap: Hashable
) -> ScoringMatrix:
    """Initialize a scoring matrix from cognate sets.

    Creates a simple identity-style matrix based on the alphabets
    observed in the cognate sets.

    Args:
        cognate_sets: The cognate sets to extract alphabets from.
        gap: Gap symbol.

    Returns:
        Initial ScoringMatrix with identity-style scoring.
    """
    # Collect all unique symbols from all sequences across all sets
    # Assuming each set has same number of sequences (e.g., pairwise)
    num_seqs = len(cognate_sets[0]) if cognate_sets else 0

    # Collect symbols per sequence position
    domains = []
    for seq_idx in range(num_seqs):
        symbols = set()
        for cog_set in cognate_sets:
            if seq_idx < len(cog_set):
                symbols.update(cog_set[seq_idx])
        domains.append(sorted([gap, *symbols]))

    # Create simple identity matrix
    return ScoringMatrix.from_sequences(
        sequences=domains, match=1.0, mismatch=-0.5, gap=gap, gap_score=-1.0
    )


def _em_learning(
    cognate_sets: list[list[Sequence[Hashable]]],
    max_iter: int,
    initial_matrix: ScoringMatrix | None,
    gap: Hashable,
    **kwargs,
) -> ScoringMatrix:
    """Learn matrix using Expectation-Maximization.

    EM algorithm:
    - E-step: Align all cognate sets with current matrix
    - M-step: Update matrix scores based on observed alignments

    Args:
        cognate_sets: Cognate sets for training.
        max_iter: Number of EM iterations.
        initial_matrix: Starting matrix (if None, creates one).
        gap: Gap symbol.
        **kwargs: Additional parameters (reserved for future use).

    Returns:
        Learned ScoringMatrix.
    """
    # Initialize matrix if not provided
    if initial_matrix is None:
        matrix = _initialize_matrix(cognate_sets, gap)
    else:
        matrix = initial_matrix.copy()

    # EM iteration
    for iteration in range(max_iter):
        # E-step: Align all cognate sets with current matrix
        alignments = []
        for cog_set in cognate_sets:
            # Get best alignment for this cognate set
            alms = multi_align(cog_set, k=1, matrix=matrix)
            if alms:
                alignments.append(alms[0])

        # M-step: Update scores based on alignments
        # Count symbol co-occurrences across all alignments
        pair_counts: Counter[tuple[Hashable, ...]] = Counter()

        for alignment in alignments:
            aln_len = len(alignment.seqs[0]) if alignment.seqs else 0
            for col_idx in range(aln_len):
                column = tuple(seq[col_idx] for seq in alignment.seqs)
                pair_counts[column] += 1

        # Update matrix scores proportional to observed frequencies
        # Use log-odds style scoring
        total_count = sum(pair_counts.values())
        if total_count > 0:
            for pair, count in pair_counts.items():
                # Simple frequency-based score
                # Normalize to get probabilities, then take log-odds
                freq = count / total_count
                score = np.log(freq + 1e-10)  # Add small constant to avoid log(0)
                matrix.scores[pair] = float(score)

    return matrix


def _gradient_descent_learning(
    cognate_sets: list[list[Sequence[Hashable]]],
    max_iter: int,
    initial_matrix: ScoringMatrix | None,
    gap: Hashable,
    **kwargs,
) -> ScoringMatrix:
    """Learn matrix using gradient descent via scipy.optimize.

    Uses L-BFGS-B optimization to find matrix that maximizes total
    alignment scores across all cognate sets.

    Args:
        cognate_sets: Cognate sets for training.
        max_iter: Maximum optimization iterations.
        initial_matrix: Starting matrix (if None, creates one).
        gap: Gap symbol.
        **kwargs: Additional parameters passed to scipy.optimize.minimize.

    Returns:
        Learned ScoringMatrix.
    """
    # Initialize matrix if not provided
    if initial_matrix is None:
        matrix = _initialize_matrix(cognate_sets, gap)
    else:
        matrix = initial_matrix.copy()

    # Get score keys in consistent order for flattening
    score_keys = sorted(matrix.scores.keys())
    num_params = len(score_keys)

    def _flatten_matrix(mat: ScoringMatrix) -> np.ndarray:
        """Flatten matrix scores to array for optimization."""
        return np.array([mat.scores[key] for key in score_keys])

    def _unflatten_to_matrix(params: np.ndarray) -> ScoringMatrix:
        """Reconstruct matrix from flattened parameters."""
        new_scores = {key: float(params[i]) for i, key in enumerate(score_keys)}
        return ScoringMatrix(
            scores=new_scores, domains=matrix.domains, gap=matrix.gap, impute_method=None
        )

    def _objective(params: np.ndarray) -> float:
        """Objective function: negative sum of alignment scores (minimize)."""
        mat = _unflatten_to_matrix(params)
        total_score = 0.0

        for cog_set in cognate_sets:
            alms = multi_align(cog_set, k=1, matrix=mat)
            if alms:
                total_score += alms[0].score

        # Return negative (we minimize, but want to maximize score)
        return -total_score

    # Initial parameters
    x0 = _flatten_matrix(matrix)

    # Optimize
    result = minimize(
        _objective,
        x0,
        method="L-BFGS-B",
        options={"maxiter": max_iter, "disp": False},
        **kwargs,
    )

    # Return optimized matrix
    return _unflatten_to_matrix(result.x)
