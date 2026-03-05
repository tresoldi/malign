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
) -> ScoringMatrix:
    """Learn a scoring matrix from cognate sets.

    Given collections of related sequences (cognates), learns a scoring matrix
    that maximizes alignment quality across all sets.

    Args:
        cognate_sets: List of cognate sets, where each set is a list of sequences.
        method: Learning method - "em" or "gradient_descent" (default: "em").
        max_iter: Maximum iterations for learning (default: 10).
        initial_matrix: Starting matrix (if None, creates identity matrix).
        gap: Gap symbol (default: "-").
        convergence_threshold: Relative score change threshold (default: 0.001).
        matrix_threshold: Frobenius norm threshold for matrix convergence (default: 0.01).
        patience: Early stopping patience (default: 5).
        bounds: Parameter bounds for gradient descent (default: (-10.0, 10.0)).
        verbose: Print convergence information (default: False).

    Returns:
        Learned ScoringMatrix optimized for the provided cognates.
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
        )
    raise ValueError(f"Unknown learning method: {method}. Use 'em' or 'gradient_descent'.")


def _initialize_matrix(
    cognate_sets: list[list[Sequence[Hashable]]],
    gap: Hashable,
) -> ScoringMatrix:
    """Initialize a scoring matrix from cognate sets.

    Args:
        cognate_sets: The cognate sets to extract alphabets from.
        gap: Gap symbol.

    Returns:
        Initial ScoringMatrix with identity-style scoring.
    """
    num_seqs = len(cognate_sets[0]) if cognate_sets else 0

    domains: list[list[Hashable]] = []
    for seq_idx in range(num_seqs):
        symbols: set[Hashable] = set()
        for cog_set in cognate_sets:
            if seq_idx < len(cog_set):
                symbols.update(cog_set[seq_idx])
        domains.append(sorted([gap, *symbols]))

    return ScoringMatrix.from_sequences(
        sequences=domains,
        match=1.0,
        mismatch=-0.5,
        gap=gap,
        gap_score=-1.0,
    )


def _sort_key(k):
    """Sort key that handles None values by treating them as empty strings."""
    if isinstance(k, tuple):
        return tuple(str(x) if x is not None else "" for x in k)
    return str(k) if k is not None else ""


def _em_learning(
    cognate_sets: list[list[Sequence[Hashable]]],
    max_iter: int,
    initial_matrix: ScoringMatrix | None,
    gap: Hashable,
    convergence_threshold: float = 0.001,
    matrix_threshold: float = 0.01,
    patience: int = 5,
    verbose: bool = False,
) -> ScoringMatrix:
    """Learn matrix using Expectation-Maximization.

    Each iteration creates a new ScoringMatrix rather than mutating.

    Args:
        cognate_sets: Cognate sets for training.
        max_iter: Maximum number of EM iterations.
        initial_matrix: Starting matrix (if None, creates one).
        gap: Gap symbol.
        convergence_threshold: Relative score change threshold (default: 0.001).
        matrix_threshold: Frobenius norm threshold (default: 0.01).
        patience: Early stopping patience (default: 5).
        verbose: Print convergence information (default: False).

    Returns:
        Learned ScoringMatrix.
    """
    matrix = _initialize_matrix(cognate_sets, gap) if initial_matrix is None else initial_matrix

    # Track convergence
    prev_total_score = None
    best_score = float("-inf")
    patience_counter = 0

    score_keys = sorted(matrix.scores.keys(), key=_sort_key)

    for iteration in range(max_iter):
        prev_matrix_values = np.array([matrix.scores[key] for key in score_keys])

        # E-step: Align all cognate sets with current matrix
        alignments = []
        total_score = 0.0

        for cog_set in cognate_sets:
            alms = align(cog_set, k=1, matrix=matrix)
            if alms:
                alignments.append(alms[0])
                if alms[0].score is not None:
                    total_score += alms[0].score

        # M-step: Build new scores from alignment co-occurrences
        pair_counts: Counter[tuple[Hashable, ...]] = Counter()

        for alignment in alignments:
            aln_len = len(alignment.seqs[0]) if alignment.seqs else 0
            for col_idx in range(aln_len):
                column = tuple(seq[col_idx] for seq in alignment.seqs)
                pair_counts[column] += 1

        # Create new matrix with updated scores
        total_count = sum(pair_counts.values())
        if total_count > 0:
            new_scores = dict(matrix.scores)
            for pair, count in pair_counts.items():
                freq = count / total_count
                score = np.log(freq + 1e-10)
                new_scores[pair] = float(score)

            matrix = ScoringMatrix(
                scores=new_scores,
                domains=matrix.domains,
                gap=matrix.gap,
                impute_method=None,
            )

        # Check convergence criteria
        converged = False

        # 1. Score-based convergence
        if prev_total_score is not None and abs(prev_total_score) > 1e-10:
            relative_change = abs(total_score - prev_total_score) / abs(prev_total_score)
            if verbose:
                print(
                    f"Iteration {iteration + 1}: score={total_score:.4f}, relative_change={relative_change:.6f}"
                )
            if relative_change < convergence_threshold:
                if verbose:
                    print(
                        f"Converged: score change {relative_change:.6f} < {convergence_threshold}"
                    )
                converged = True
        elif verbose:
            print(f"Iteration {iteration + 1}: score={total_score:.4f}")

        # 2. Matrix-based convergence
        current_matrix_values = np.array([matrix.scores[key] for key in score_keys])
        frobenius_norm = float(np.linalg.norm(current_matrix_values - prev_matrix_values))
        if verbose:
            print(f"  Matrix Frobenius norm: {frobenius_norm:.6f}")
        if frobenius_norm < matrix_threshold:
            if verbose:
                print(f"Converged: Frobenius norm {frobenius_norm:.6f} < {matrix_threshold}")
            converged = True

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
) -> ScoringMatrix:
    """Learn matrix using gradient descent via scipy.optimize.

    Args:
        cognate_sets: Cognate sets for training.
        max_iter: Maximum optimization iterations.
        initial_matrix: Starting matrix (if None, creates one).
        gap: Gap symbol.
        bounds: Parameter bounds as (min, max) tuple (default: (-10.0, 10.0)).
        patience: Early stopping patience (default: 5).
        verbose: Print optimization information (default: False).

    Returns:
        Learned ScoringMatrix.
    """
    matrix = _initialize_matrix(cognate_sets, gap) if initial_matrix is None else initial_matrix

    score_keys = sorted(matrix.scores.keys(), key=_sort_key)
    num_params = len(score_keys)

    # Early stopping tracking
    best_objective = float("inf")
    patience_counter = 0
    iteration_count = [0]

    def _flatten_matrix(mat: ScoringMatrix) -> np.ndarray:
        return np.array([mat.scores[key] for key in score_keys])

    def _unflatten_to_matrix(params: np.ndarray) -> ScoringMatrix:
        new_scores = {key: float(params[i]) for i, key in enumerate(score_keys)}
        return ScoringMatrix(
            scores=new_scores,
            domains=matrix.domains,
            gap=matrix.gap,
            impute_method=None,
        )

    def _objective(params: np.ndarray) -> float:
        mat = _unflatten_to_matrix(params)
        total_score = 0.0

        for cog_set in cognate_sets:
            alms = align(cog_set, k=1, matrix=mat)
            if alms and alms[0].score is not None:
                total_score += alms[0].score

        return -total_score

    def _callback(params: np.ndarray) -> bool:
        nonlocal best_objective, patience_counter

        iteration_count[0] += 1
        current_objective = _objective(params)

        if verbose:
            print(f"Iteration {iteration_count[0]}: objective={current_objective:.4f}")

        if current_objective < best_objective:
            best_objective = current_objective
            patience_counter = 0
            if verbose:
                print(f"  New best objective: {best_objective:.4f}")
        else:
            patience_counter += 1
            if verbose:
                print(f"  No improvement: patience {patience_counter}/{patience}")

            if patience_counter >= patience:
                if verbose:
                    print(f"Early stopping: no improvement for {patience} iterations")
                return True

        return False

    x0 = _flatten_matrix(matrix)
    param_bounds = [bounds] * num_params

    result = minimize(
        _objective,
        x0,
        method="L-BFGS-B",
        bounds=param_bounds,
        callback=_callback,
        options={"maxiter": max_iter, "disp": verbose},
    )

    if verbose:
        print(f"Optimization finished: {result.message}")
        print(f"Final objective: {result.fun:.4f}")
        print(f"Total iterations: {iteration_count[0]}")

    return _unflatten_to_matrix(result.x)
