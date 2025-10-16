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

from .malign import align
from .scoring_matrix import ScoringMatrix


def learn_matrix(
    cognate_sets: list[list[Sequence[Hashable]]],
    method: str = "em",
    max_iter: int = 10,
    initial_matrix: ScoringMatrix | None = None,
    gap: Hashable = "-",
    convergence_threshold: float = 0.001,
    matrix_threshold: float = 0.01,
    patience: int = 5,
    bounds: tuple[float, float] = (-10.0, 10.0),
    verbose: bool = False,
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
        convergence_threshold: Relative score change threshold for convergence (default: 0.001).
            Stops if |new_score - old_score| / |old_score| < threshold.
        matrix_threshold: Frobenius norm threshold for matrix convergence (default: 0.01).
            Stops if sqrt(sum((new[k] - old[k])**2)) < threshold.
        patience: Early stopping patience - stop after N iterations without improvement (default: 5).
        bounds: Parameter bounds for gradient descent as (min, max) (default: (-10.0, 10.0)).
        verbose: Print convergence information during learning (default: False).
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
        Phase 3.7 adds convergence detection and early stopping for both EM and gradient descent.
    """
    if method == "em":
        return _em_learning(
            cognate_sets,
            max_iter,
            initial_matrix,
            gap,
            convergence_threshold=convergence_threshold,
            matrix_threshold=matrix_threshold,
            patience=patience,
            verbose=verbose,
            **kwargs,
        )
    if method == "gradient_descent":
        return _gradient_descent_learning(
            cognate_sets,
            max_iter,
            initial_matrix,
            gap,
            bounds=bounds,
            patience=patience,
            verbose=verbose,
            **kwargs,
        )
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
    convergence_threshold: float = 0.001,
    matrix_threshold: float = 0.01,
    patience: int = 5,
    verbose: bool = False,
    **kwargs,
) -> ScoringMatrix:
    """Learn matrix using Expectation-Maximization.

    EM algorithm:
    - E-step: Align all cognate sets with current matrix
    - M-step: Update matrix scores based on observed alignments

    Phase 3.7 enhancements:
    - Convergence detection (score-based and matrix-based)
    - Early stopping with patience

    Args:
        cognate_sets: Cognate sets for training.
        max_iter: Maximum number of EM iterations.
        initial_matrix: Starting matrix (if None, creates one).
        gap: Gap symbol.
        convergence_threshold: Relative score change threshold (default: 0.001).
        matrix_threshold: Frobenius norm threshold (default: 0.01).
        patience: Early stopping patience (default: 5).
        verbose: Print convergence information (default: False).
        **kwargs: Additional parameters (reserved for future use).

    Returns:
        Learned ScoringMatrix.
    """
    # Initialize matrix if not provided
    if initial_matrix is None:
        matrix = _initialize_matrix(cognate_sets, gap)
    else:
        matrix = initial_matrix.copy()

    # Track convergence
    prev_total_score = None
    best_score = float("-inf")
    patience_counter = 0

    # Get score keys for matrix comparison
    # Use custom key to handle None values in tuples
    def _sort_key(k):
        """Sort key that handles None values by treating them as empty strings."""
        if isinstance(k, tuple):
            return tuple(str(x) if x is not None else "" for x in k)
        return str(k) if k is not None else ""

    score_keys = sorted(matrix.scores.keys(), key=_sort_key)

    # EM iteration
    for iteration in range(max_iter):
        # Store previous matrix for convergence check
        prev_matrix_values = np.array([matrix.scores[key] for key in score_keys])

        # E-step: Align all cognate sets with current matrix
        alignments = []
        total_score = 0.0

        for cog_set in cognate_sets:
            # Get best alignment for this cognate set
            alms = align(cog_set, k=1, matrix=matrix)
            if alms:
                alignments.append(alms[0])
                if alms[0].score is not None:
                    total_score += alms[0].score

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

        # Check convergence criteria
        converged = False

        # 1. Score-based convergence: relative change in total score
        if prev_total_score is not None and abs(prev_total_score) > 1e-10:
            relative_change = abs(total_score - prev_total_score) / abs(prev_total_score)
            if verbose:
                print(f"Iteration {iteration + 1}: score={total_score:.4f}, relative_change={relative_change:.6f}")
            if relative_change < convergence_threshold:
                if verbose:
                    print(f"Converged: score change {relative_change:.6f} < {convergence_threshold}")
                converged = True
        elif verbose:
            print(f"Iteration {iteration + 1}: score={total_score:.4f}")

        # 2. Matrix-based convergence: Frobenius norm of parameter change
        current_matrix_values = np.array([matrix.scores[key] for key in score_keys])
        frobenius_norm = float(np.linalg.norm(current_matrix_values - prev_matrix_values))
        if verbose:
            print(f"  Matrix Frobenius norm: {frobenius_norm:.6f}")
        if frobenius_norm < matrix_threshold:
            if verbose:
                print(f"Converged: Frobenius norm {frobenius_norm:.6f} < {matrix_threshold}")
            converged = True

        # Stop if converged (OR logic - either criterion met)
        if converged:
            if verbose:
                print(f"Stopping early at iteration {iteration + 1}")
            break

        # 3. Early stopping with patience
        if total_score > best_score:
            best_score = total_score
            patience_counter = 0
        else:
            patience_counter += 1
            if verbose:
                print(f"  No improvement: patience {patience_counter}/{patience}")
            if patience_counter >= patience:
                if verbose:
                    print(f"Early stopping: no improvement for {patience} iterations")
                break

        # Update previous score for next iteration
        prev_total_score = total_score

    return matrix


def _gradient_descent_learning(
    cognate_sets: list[list[Sequence[Hashable]]],
    max_iter: int,
    initial_matrix: ScoringMatrix | None,
    gap: Hashable,
    bounds: tuple[float, float] = (-10.0, 10.0),
    patience: int = 5,
    verbose: bool = False,
    **kwargs,
) -> ScoringMatrix:
    """Learn matrix using gradient descent via scipy.optimize.

    Uses L-BFGS-B optimization to find matrix that maximizes total
    alignment scores across all cognate sets.

    Phase 3.7 enhancements:
    - Parameter bounds to prevent extreme values
    - Early stopping with patience

    Args:
        cognate_sets: Cognate sets for training.
        max_iter: Maximum optimization iterations.
        initial_matrix: Starting matrix (if None, creates one).
        gap: Gap symbol.
        bounds: Parameter bounds as (min, max) tuple (default: (-10.0, 10.0)).
        patience: Early stopping patience (default: 5).
        verbose: Print optimization information (default: False).
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
    # Use custom key to handle None values in tuples
    def _sort_key(k):
        """Sort key that handles None values by treating them as empty strings."""
        if isinstance(k, tuple):
            return tuple(str(x) if x is not None else "" for x in k)
        return str(k) if k is not None else ""

    score_keys = sorted(matrix.scores.keys(), key=_sort_key)
    num_params = len(score_keys)

    # Early stopping tracking
    best_objective = float("inf")
    patience_counter = 0
    iteration_count = [0]  # Use list to allow modification in nested function

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
            alms = align(cog_set, k=1, matrix=mat)
            if alms and alms[0].score is not None:
                total_score += alms[0].score

        # Return negative (we minimize, but want to maximize score)
        return -total_score

    def _callback(params: np.ndarray) -> bool:
        """Callback for early stopping with patience.

        Returns:
            True to stop optimization, False to continue.
        """
        nonlocal best_objective, patience_counter

        iteration_count[0] += 1
        current_objective = _objective(params)

        if verbose:
            print(f"Iteration {iteration_count[0]}: objective={current_objective:.4f}")

        # Check if this is the best objective so far
        if current_objective < best_objective:
            best_objective = current_objective
            patience_counter = 0
            if verbose:
                print(f"  New best objective: {best_objective:.4f}")
        else:
            patience_counter += 1
            if verbose:
                print(f"  No improvement: patience {patience_counter}/{patience}")

            # Stop if patience exceeded
            if patience_counter >= patience:
                if verbose:
                    print(f"Early stopping: no improvement for {patience} iterations")
                return True  # Stop optimization

        return False  # Continue optimization

    # Initial parameters
    x0 = _flatten_matrix(matrix)

    # Set bounds for all parameters
    param_bounds = [bounds] * num_params

    # Optimize with bounds and callback
    result = minimize(
        _objective,
        x0,
        method="L-BFGS-B",
        bounds=param_bounds,
        callback=_callback,
        options={"maxiter": max_iter, "disp": verbose},
        **kwargs,
    )

    if verbose:
        print(f"Optimization finished: {result.message}")
        print(f"Final objective: {result.fun:.4f}")
        print(f"Total iterations: {iteration_count[0]}")

    # Return optimized matrix
    return _unflatten_to_matrix(result.x)
