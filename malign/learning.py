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

from .blocks import _adjust_pair_counts_for_blocks
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
    prior_matrix: ScoringMatrix | None = None,
    prior_weight: float = 0.5,
) -> ScoringMatrix:
    """Learn a scoring matrix from cognate sets.

    Given collections of related sequences (cognates), learns a scoring matrix
    that maximizes alignment quality across all sets.

    Args:
        cognate_sets: List of cognate sets, where each set is a list of sequences.
        method: Learning method - "em" or "gradient_descent" (default: "em").
        max_iter: Maximum iterations for learning (default: 10).
        initial_matrix: Starting matrix (if None, creates identity matrix;
            if prior_matrix is given, uses it as initial).
        gap: Gap symbol (default: "-").
        convergence_threshold: Relative score change threshold (default: 0.001).
        matrix_threshold: Frobenius norm threshold for matrix convergence (default: 0.01).
        patience: Early stopping patience (default: 5).
        bounds: Parameter bounds for gradient descent (default: (-10.0, 10.0)).
        verbose: Print convergence information (default: False).
        prior_matrix: Optional prior matrix (e.g. from ``ScoringMatrix.from_distfeat()``).
            Blended with data-driven scores during EM learning. Symbols in the prior
            but absent from the cognate sets are included in the output.
        prior_weight: Initial regularization strength for the prior (default: 0.5).
            Decays linearly to 0 over ``max_iter`` iterations (EM only).

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
            prior_matrix=prior_matrix,
            prior_weight=prior_weight,
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
    prior_matrix: ScoringMatrix | None = None,
    prior_weight: float = 0.5,
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
        prior_matrix: Optional prior matrix for regularization.
        prior_weight: Initial prior blending weight (decays linearly).

    Returns:
        Learned ScoringMatrix.
    """
    if initial_matrix is not None:
        matrix = initial_matrix
    elif prior_matrix is not None:
        matrix = prior_matrix
    else:
        matrix = _initialize_matrix(cognate_sets, gap)

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

            # Prior regularization with linear decay
            alpha = prior_weight * (1.0 - iteration / max_iter)
            if prior_matrix is not None and alpha > 0:
                for pair_key in new_scores:
                    prior_score = prior_matrix.scores.get(pair_key)
                    if prior_score is not None:
                        new_scores[pair_key] += alpha * prior_score

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


def bootstrap_matrix(
    pairs: list[tuple[Sequence[Hashable], Sequence[Hashable]]],
    max_iter: int = 20,
    initial_matrix: ScoringMatrix | None = None,
    gap: Hashable = "-",
    alignment_method: str = "anw",
    convergence_threshold: float = 0.001,
    matrix_threshold: float = 0.01,
    patience: int = 5,
    gap_score: float = -1.0,
    verbose: bool = False,
    prior_matrix: ScoringMatrix | None = None,
    prior_weight: float = 0.5,
    block_merge: bool = False,
    max_block_size: int = 2,
) -> ScoringMatrix:
    """Learn a scoring matrix from sequence pairs without pre-clustered cognate sets.

    Iteratively aligns pairs with the current matrix, counts column
    co-occurrences, and re-estimates scores using log-odds. Unlike
    ``learn_matrix()`` (which requires pre-clustered cognate sets), this
    function works with arbitrary sequence pairs.

    When ``prior_matrix`` is provided, each M-step blends the data-driven
    log-odds scores with the prior using additive regularization.  The
    blending weight (``prior_weight``) decays linearly to zero over
    ``max_iter`` iterations, so the prior dominates early (when data signal
    is weak) and fades as alignments stabilize.

    When ``block_merge`` is True, complementary-gap block columns are
    detected in each alignment and the gap-containing columns within
    blocks are decremented from pair counts, reducing gap inflation
    caused by diphthongization and metathesis patterns.

    Args:
        pairs: List of (sequence_a, sequence_b) tuples.
        max_iter: Maximum number of bootstrap iterations (default: 20).
        initial_matrix: Starting matrix. If None, creates an identity matrix.
        gap: Gap symbol (default: "-").
        alignment_method: Alignment method for each iteration (default: "anw").
        convergence_threshold: Relative score change threshold (default: 0.001).
        matrix_threshold: Frobenius norm threshold for convergence (default: 0.01).
        patience: Early stopping patience (default: 5).
        gap_score: Score for gap alignment in initial matrix (default: -1.0).
        verbose: Print convergence information (default: False).
        prior_matrix: Optional phonological/feature prior (e.g. from
            ``ScoringMatrix.from_distfeat()``). Symbols in the prior but
            absent from the pairs are included in the output matrix.
        prior_weight: Initial regularization strength for the prior
            (default: 0.5). Decays linearly to 0 over ``max_iter``.
        block_merge: If True, adjust pair counts by detecting
            complementary-gap blocks and reducing gap inflation (default: False).
        max_block_size: Maximum block size for block detection (default: 2).
            Only used when block_merge is True.

    Returns:
        Learned ScoringMatrix.

    Raises:
        ValueError: If pairs is empty.
    """
    if not pairs:
        raise ValueError("At least one pair is required for bootstrap learning.")

    # Collect all symbols from both sides
    left_symbols: set[Hashable] = set()
    right_symbols: set[Hashable] = set()
    for seq_a, seq_b in pairs:
        left_symbols.update(seq_a)
        right_symbols.update(seq_b)

    # Expand symbol sets with prior matrix domains
    if prior_matrix is not None:
        if len(prior_matrix.domains) >= 1:
            left_symbols |= set(prior_matrix.domains[0])
        if len(prior_matrix.domains) >= 2:
            right_symbols |= set(prior_matrix.domains[1])
        left_symbols.discard(gap)
        right_symbols.discard(gap)

    # Initialize matrix
    if initial_matrix is not None:
        matrix = initial_matrix
    elif prior_matrix is not None:
        matrix = prior_matrix
    else:
        matrix = ScoringMatrix.from_sequences(
            sequences=[sorted([gap, *left_symbols]), sorted([gap, *right_symbols])],
            match=1.0,
            mismatch=-0.5,
            gap=gap,
            gap_score=gap_score,
        )

    score_keys = sorted(matrix.scores.keys(), key=_sort_key)
    prev_total_score = None
    best_score = float("-inf")
    patience_counter = 0

    for iteration in range(max_iter):
        prev_matrix_values = np.array([matrix.scores[key] for key in score_keys])

        # E-step: Align all pairs with current matrix
        pair_counts: Counter[tuple[Hashable, ...]] = Counter()
        total_score = 0.0
        iteration_alignments: list = []

        for seq_a, seq_b in pairs:
            alms = align(
                [list(seq_a), list(seq_b)],
                method=alignment_method,
                k=1,
                matrix=matrix,
            )

            if alms and alms[0].score is not None:
                total_score += alms[0].score
                iteration_alignments.append(alms[0])
                aln_len = len(alms[0].seqs[0])
                for col_idx in range(aln_len):
                    column = tuple(seq[col_idx] for seq in alms[0].seqs)
                    pair_counts[column] += 1

        # Block adjustment: reduce gap inflation from complementary-gap patterns
        if block_merge:
            for aln in iteration_alignments:
                _adjust_pair_counts_for_blocks(pair_counts, aln, gap, max_block_size)

        # M-step: Build log-odds scores with Laplace smoothing
        total_count = sum(pair_counts.values())
        if total_count > 0:
            # Marginal frequencies per domain
            all_symbols_0 = sorted(left_symbols | {gap})
            all_symbols_1 = sorted(right_symbols | {gap})

            marginal_0: Counter[Hashable] = Counter()
            marginal_1: Counter[Hashable] = Counter()
            for (s0, s1), count in pair_counts.items():
                marginal_0[s0] += count
                marginal_1[s1] += count

            num_symbol_pairs = len(all_symbols_0) * len(all_symbols_1)
            smoothed_total = total_count + num_symbol_pairs  # Laplace

            new_scores: dict[tuple[Hashable, ...], float] = {}
            for s0 in all_symbols_0:
                for s1 in all_symbols_1:
                    pair_key = (s0, s1)
                    # Observed frequency with Laplace smoothing (+1)
                    observed = (pair_counts.get(pair_key, 0) + 1) / smoothed_total
                    # Expected frequency (independence model)
                    freq_0 = (marginal_0.get(s0, 0) + 1) / (smoothed_total)
                    freq_1 = (marginal_1.get(s1, 0) + 1) / (smoothed_total)
                    expected = freq_0 * freq_1

                    if expected > 0:
                        new_scores[pair_key] = float(np.log(observed / expected))
                    else:
                        new_scores[pair_key] = 0.0

            # Gap-involving scores must be non-positive: gaps are penalties,
            # never rewards.  Without this clamp the log-odds ratio can go
            # positive for rare symbols, creating a feedback loop where the
            # aligner inserts gratuitous gaps that inflate alignment length.
            for pair_key in new_scores:
                if gap in pair_key:
                    new_scores[pair_key] = min(new_scores[pair_key], 0.0)

            # All-gap vector always scores 0
            new_scores[(gap, gap)] = 0.0

            # Prior regularization with linear decay
            alpha = prior_weight * (1.0 - iteration / max_iter)
            if prior_matrix is not None and alpha > 0:
                for pair_key in new_scores:
                    prior_score = prior_matrix.scores.get(pair_key)
                    if prior_score is not None:
                        new_scores[pair_key] += alpha * prior_score

            matrix = ScoringMatrix(
                scores=new_scores,
                domains=[all_symbols_0, all_symbols_1],
                gap=gap,
                impute_method=None,
            )

            # Update score_keys for convergence tracking
            score_keys = sorted(matrix.scores.keys(), key=_sort_key)

        # Check convergence
        converged = False

        if prev_total_score is not None and abs(prev_total_score) > 1e-10:
            relative_change = abs(total_score - prev_total_score) / abs(prev_total_score)
            if verbose:
                print(
                    f"Iteration {iteration + 1}: score={total_score:.4f}, "
                    f"relative_change={relative_change:.6f}"
                )
            if relative_change < convergence_threshold:
                if verbose:
                    print(
                        f"Converged: score change {relative_change:.6f} < {convergence_threshold}"
                    )
                converged = True
        elif verbose:
            print(f"Iteration {iteration + 1}: score={total_score:.4f}")

        # Matrix-based convergence
        current_matrix_values = np.array([matrix.scores.get(key, 0.0) for key in score_keys])
        if len(prev_matrix_values) == len(current_matrix_values):
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

        # Early stopping with patience
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
