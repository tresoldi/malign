"""Module for scoring matrices."""

import itertools
import math
from collections import Counter
from collections.abc import Hashable, Iterable, Sequence
from dataclasses import dataclass
from typing import cast

import numpy as np
import yaml
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.experimental import enable_iterative_imputer  # noqa: F401
from sklearn.impute import IterativeImputer, SimpleImputer
from sklearn.linear_model import BayesianRidge
from sklearn.neighbors import KNeighborsRegressor
from sklearn.tree import DecisionTreeRegressor
from tabulate import tabulate


def _symbol_sort_key(symbol: Hashable) -> tuple[str, str]:
    """Deterministic sort key that supports heterogeneous hashable symbols."""
    return (type(symbol).__name__, repr(symbol))


def _build_domains(
    domains: list[list[Hashable]] | tuple[tuple[Hashable, ...], ...] | None,
    scores: dict[tuple[Hashable, ...], float],
    gap: Hashable,
) -> tuple[tuple[Hashable, ...], ...]:
    """Build sorted, deduplicated domains from user input or scores.

    Args:
        domains: User-provided domains or None to infer from scores.
        scores: Score dictionary to infer domains from if not provided.
        gap: Gap symbol to ensure is in each domain.

    Returns:
        Tuple of tuples of domain symbols.
    """
    raw: Sequence[Iterable[Hashable]]
    if domains is not None:
        raw = domains
    else:
        inferred: list[set[Hashable]] = [
            {symbol for symbol in domain if symbol is not None} | {gap}
            for domain in zip(*scores.keys(), strict=False)
        ]
        raw = inferred

    return tuple(tuple(sorted(set(d), key=_symbol_sort_key)) for d in raw)


def _validate_scores(
    scores: dict[tuple[Hashable, ...], float],
    domains: tuple[tuple[Hashable, ...], ...],
) -> None:
    """Validate that score keys are consistent and use domain symbols.

    Args:
        scores: Score dictionary to validate.
        domains: Domain tuples to validate against.

    Raises:
        ValueError: If scores have inconsistent key lengths or unknown symbols.
    """
    key_lens = {len(key) for key in scores}
    if len(key_lens) > 1:
        raise ValueError("Different domain-lengths in `scores`.")

    scores_domains = zip(*scores.keys(), strict=False)
    found = [
        all(symbol in ref_domain for symbol in scores_domain if symbol is not None)
        for scores_domain, ref_domain in zip(scores_domains, domains, strict=False)
    ]
    if not all(found):
        raise ValueError("`scores` has symbols not in domains.")


def _fill_matrix(
    scores: dict[tuple[Hashable, ...], float],
    domains: tuple[tuple[Hashable, ...], ...],
    impute_method: str,
) -> dict[tuple[Hashable, ...], float]:
    """Fill a sparse matrix using imputation, returning a new scores dict.

    Args:
        scores: Existing scores to train the imputer on.
        domains: Domain tuples.
        impute_method: Imputation method name.

    Returns:
        New scores dict with imputed values filled in.
    """
    encoder = list(
        itertools.chain.from_iterable(
            [(domain_idx, value) for value in [None, *domain_values]]
            for domain_idx, domain_values in enumerate(domains)
        )
    )

    domains_with_none = [[None, *domain] for domain in domains]
    train_matrix = []
    imp_matrix = []
    for cat_vector in itertools.product(*domains_with_none):
        if len([val for val in cat_vector if val]) < 2:
            continue

        mh_vector = [cat_vector[element[0]] == element[1] for element in encoder]
        score = scores.get(cat_vector)
        if score is not None:
            train_matrix.append([*mh_vector, score])
        else:
            imp_matrix.append([*mh_vector, np.nan])

    if not imp_matrix:
        return scores

    if impute_method == "decision_tree":
        estimator = DecisionTreeRegressor(max_features="sqrt", random_state=0)
        imputer = IterativeImputer(missing_values=np.nan, random_state=0, estimator=estimator)
    elif impute_method == "extra_trees":
        estimator = ExtraTreesRegressor(n_estimators=10, random_state=0)
        imputer = IterativeImputer(missing_values=np.nan, random_state=0, estimator=estimator)
    elif impute_method == "k_neighbors":
        n_neighbors = min(5, len(train_matrix) - 1) if len(train_matrix) > 1 else 1
        estimator = KNeighborsRegressor(n_neighbors=n_neighbors)
        imputer = IterativeImputer(missing_values=np.nan, random_state=0, estimator=estimator)
    elif impute_method == "bayesian_ridge":
        estimator = BayesianRidge()
        imputer = IterativeImputer(missing_values=np.nan, random_state=0, estimator=estimator)
    elif impute_method in ["mean", "median"]:
        imputer = SimpleImputer(missing_values=np.nan, strategy=impute_method)
    else:
        imputer = SimpleImputer(missing_values=np.nan, strategy="mean")

    imputer.fit(train_matrix)
    trans_matrix = imputer.transform(imp_matrix)

    new_scores = dict(scores)
    for row in trans_matrix:
        mh_vector, value = row[:-1], row[-1]
        cat_vector = tuple(encoder[idx][1] for idx, flag in enumerate(mh_vector) if flag)
        new_scores[cat_vector] = float(value)

    return new_scores


@dataclass(frozen=True)
class ScoringMatrix:
    """Frozen dataclass for sequence alignment scoring matrices.

    A scoring matrix stores alignment scores for symbol tuples across multiple
    domains. It supports asymmetric scoring and treats gaps as normal symbols.

    Use factory classmethods to create instances:
    - ``ScoringMatrix(scores=..., ...)`` for direct construction
    - ``ScoringMatrix.from_yaml(filename)`` to load from YAML
    - ``ScoringMatrix.from_sequences(sequences)`` for identity-style matrices
    """

    scores: dict[tuple[Hashable, ...], float]
    domains: tuple[tuple[Hashable, ...], ...]
    gap: Hashable = "-"
    num_domains: int = 0

    def __init__(
        self,
        scores: dict[tuple[Hashable, ...], float] | None = None,
        domains: list[list[Hashable]] | tuple[tuple[Hashable, ...], ...] | None = None,
        gap: Hashable = "-",
        impute_method: str | None = "mean",
    ):
        """Initialize a scoring matrix.

        Args:
            scores: Scoring dictionary mapping symbol tuples to float scores.
            domains: Domain symbol lists. Inferred from scores if not provided.
            gap: Gap symbol. Defaults to "-".
            impute_method: Method for filling sparse matrices. Choices: "mean",
                "median", "decision_tree", "extra_trees", "k_neighbors",
                "bayesian_ridge". None disables imputation. Defaults to "mean".
        """
        object.__setattr__(self, "gap", gap)

        if scores is None:
            # Empty matrix - will be populated via from_yaml
            object.__setattr__(self, "scores", {})
            object.__setattr__(self, "domains", ())
            object.__setattr__(self, "num_domains", 0)
            return

        built_domains = _build_domains(domains, scores, gap)
        _validate_scores(scores, built_domains)

        num = len(built_domains)
        new_scores = dict(scores)
        new_scores[tuple([gap] * num)] = 0.0

        if impute_method:
            if impute_method not in [
                "mean",
                "median",
                "decision_tree",
                "extra_trees",
                "k_neighbors",
                "bayesian_ridge",
            ]:
                raise ValueError(f"Unknown imputation method: {impute_method}.")
            new_scores = _fill_matrix(new_scores, built_domains, impute_method)

        object.__setattr__(self, "scores", new_scores)
        object.__setattr__(self, "domains", built_domains)
        object.__setattr__(self, "num_domains", num)

    @classmethod
    def from_yaml(cls, filename: str, impute_method: str | None = "mean") -> "ScoringMatrix":
        """Load a scoring matrix from a YAML file.

        Args:
            filename: Path to the YAML file.
            impute_method: Imputation method for filling sparse matrices.

        Returns:
            A new ScoringMatrix loaded from the file.
        """
        with open(filename, encoding="utf-8") as yaml_handler:
            serial_data = yaml.safe_load(yaml_handler)

        gap = serial_data["gap"]
        domains = [list(d) for d in serial_data["domains"]]
        scores = {
            tuple(None if sub_key == "NULL" else sub_key for sub_key in key.split(" / ")): float(
                value
            )
            for key, value in serial_data["scores"].items()
        }

        return cls(scores=scores, domains=domains, gap=gap, impute_method=impute_method)

    @classmethod
    def from_sequences(
        cls,
        sequences: list[list[Hashable]],
        match: float = 1.0,
        mismatch: float = -1.0,
        gap: Hashable = "-",
        gap_score: float = -1.0,
        impute_method: str | None = "mean",
    ) -> "ScoringMatrix":
        """Create a scoring matrix from sequences with match/mismatch scoring.

        Args:
            sequences: List of sequence alphabets.
            match: Score for matching symbols (default: 1.0).
            mismatch: Score for mismatching symbols (default: -1.0).
            gap: Gap symbol (default: "-").
            gap_score: Score for gaps (default: -1.0).
            impute_method: Method for filling missing scores (default: "mean").

        Returns:
            A new ScoringMatrix with basic match/mismatch scoring.
        """
        domains = [sorted({gap, *seq}, key=_symbol_sort_key) for seq in sequences]

        scores = {}
        for combo in itertools.product(*domains):
            if all(s == gap for s in combo):
                scores[combo] = 0.0
            elif gap in combo:
                scores[combo] = gap_score
            elif len(set(combo)) == 1:
                scores[combo] = match
            else:
                scores[combo] = mismatch

        return cls(scores=scores, domains=domains, gap=gap, impute_method=impute_method)

    @classmethod
    def from_distfeat(
        cls,
        sequences: list[list[str]],
        gap: str = "-",
        gap_score: float = -1.0,
        system: str | None = None,
        impute_method: str | None = "mean",
    ) -> "ScoringMatrix":
        """Create a scoring matrix from phonological feature distances.

        Uses the distfeat library to compute pairwise segment distances
        based on distinctive phonological features, then converts them
        to similarity scores.

        Args:
            sequences: List of sequence alphabets (lists of IPA segments).
            gap: Gap symbol (default: "-").
            gap_score: Score for gap alignments (default: -1.0).
            system: Distfeat feature system name (default: None for default system).
            impute_method: Method for filling missing scores (default: "mean").

        Returns:
            A new ScoringMatrix with feature-distance-based scores.

        Raises:
            ImportError: If distfeat is not installed.
        """
        try:
            from distfeat import distance as _distfeat_distance
        except ImportError:
            raise ImportError(
                "distfeat is required for from_distfeat(). "
                "Install it with: pip install malign[features]"
            ) from None

        domains = [sorted({gap, *seq}, key=_symbol_sort_key) for seq in sequences]

        scores: dict[tuple[Hashable, ...], float] = {}
        for combo in itertools.product(*domains):
            if all(s == gap for s in combo):
                scores[combo] = 0.0
            elif gap in combo:
                scores[combo] = gap_score
            else:
                # Convert distance to similarity: score = 1 - distance
                # distfeat.distance returns values in [0, 1]
                dist = _distfeat_distance(combo[0], combo[1], system=system)
                scores[combo] = 1.0 - dist

        hashable_domains = cast("list[list[Hashable]]", domains)
        return cls(scores=scores, domains=hashable_domains, gap=gap, impute_method=impute_method)

    @classmethod
    def from_substitution_counts(
        cls,
        counts: dict[tuple[Hashable, ...], int],
        gap: Hashable = "-",
        gap_score: float = -1.0,
        impute_method: str | None = "mean",
    ) -> "ScoringMatrix":
        """Create an asymmetric matrix from observed substitution frequencies.

        Converts raw substitution counts to log-odds scores:
        ``score = log(observed_freq / expected_freq)`` where expected
        frequency is the product of marginal frequencies (independence model).

        This is the standard approach used by BLOSUM/PAM matrices and is
        directly applicable to historical sound change rates in linguistics.

        Args:
            counts: Mapping from symbol tuples to observed counts.
                E.g. ``{("p", "b"): 15, ("b", "p"): 3}``.
            gap: Gap symbol (default: "-").
            gap_score: Score for gap alignments (default: -1.0).
            impute_method: Method for filling missing scores (default: "mean").

        Returns:
            A new ScoringMatrix with asymmetric log-odds scores.

        Raises:
            ValueError: If counts is empty or has inconsistent key lengths.
        """
        if not counts:
            raise ValueError("counts must be non-empty.")

        # Determine arity from keys
        arity = len(next(iter(counts)))

        # Collect all symbols per position
        symbols_per_pos: list[set[Hashable]] = [set() for _ in range(arity)]
        for key in counts:
            for d, sym in enumerate(key):
                symbols_per_pos[d].add(sym)

        domains = [sorted({gap, *syms}, key=_symbol_sort_key) for syms in symbols_per_pos]

        # Compute total and marginal frequencies
        total = sum(counts.values())
        marginals: list[Counter] = [Counter() for _ in range(arity)]
        for key, count in counts.items():
            for d, sym in enumerate(key):
                marginals[d][sym] += count

        # Convert to log-odds scores
        scores: dict[tuple[Hashable, ...], float] = {}
        for combo in itertools.product(*domains):
            if all(s == gap for s in combo):
                scores[combo] = 0.0
                continue
            if gap in combo:
                scores[combo] = gap_score
                continue

            observed = counts.get(combo, 0)
            if observed == 0:
                # No observation: leave for imputation
                continue

            observed_freq = observed / total
            expected_freq = math.prod(marginals[d][sym] / total for d, sym in enumerate(combo))

            if expected_freq > 0:
                scores[combo] = math.log(observed_freq / expected_freq)
            else:
                continue  # Leave for imputation

        return cls(scores=scores, domains=domains, gap=gap, impute_method=impute_method)

    def save(self, filename: str) -> None:
        """Serialize the matrix to YAML format.

        Args:
            filename: Path to the YAML output file.
        """
        if any(" / " in str(symbol) for domain in self.domains for symbol in domain):
            raise ValueError("At least one domain uses reserved symbols.")

        def _allow_none_key(obj):
            return tuple("0000000000" if k is None else k for k in obj)

        _scores = {
            " / ".join("NULL" if k is None else str(k) for k in key): float(self.scores[key])
            for key in sorted(self.scores, key=_allow_none_key)
        }

        with open(filename, "w", encoding="utf-8") as yaml_handler:
            serial_data = {
                "gap": self.gap,
                "domains": [list(d) for d in self.domains],
                "domain_range": list(range(self.num_domains)),
                "scores": _scores,
            }
            yaml.dump(serial_data, yaml_handler, default_flow_style=False, allow_unicode=True)

    def compute_submatrices(self, domains: list[tuple[int, ...]]) -> dict:
        """Compute submatrices from a collection of domain indices.

        Args:
            domains: List of domain index tuples for submatrices.

        Returns:
            Dictionary mapping domain tuples to ScoringMatrix instances.
        """
        domain_range = tuple(range(self.num_domains))
        sub_matrix_scores: dict[tuple[Hashable, ...], ScoringMatrix] = {}

        for sub_domain in domains:
            sub_scores: dict[tuple[Hashable, ...], float] = {}
            for key, value in self.scores.items():
                check = [key[idx] is not None for idx in sub_domain]
                check += [key[idx] is None for idx in domain_range if idx not in sub_domain]

                if all(check):
                    sub_key = tuple(k for k in key if k is not None)
                    sub_scores[sub_key] = value

            sub_matrix_scores[sub_domain] = ScoringMatrix(sub_scores, gap=self.gap)

        return sub_matrix_scores

    def tabulate(self) -> str:
        """Build a tabulated string representation of the matrix.

        Returns:
            A string with a tabular representation for human inspection.
        """
        rows = []
        if self.num_domains == 2:
            for symbol_a in self.domains[0]:
                row = [symbol_a] + [self.scores[symbol_a, symbol_b] for symbol_b in self.domains[1]]
                rows.append(row)
            headers = ["", *(str(symbol) for symbol in self.domains[1])]

        elif self.num_domains == 3:
            for symbol_a in self.domains[0]:
                row = [symbol_a] + [
                    self.scores.get((symbol_a, symbol_b, symbol_c), "-")
                    for symbol_b, symbol_c in itertools.product(self.domains[1], self.domains[2])
                ]
                rows.append(row)
            headers = [""] + [
                "/" + "/".join(str(symbol) for symbol in sub_key)
                for sub_key in itertools.product(self.domains[1], self.domains[2])
            ]
        else:
            raise ValueError("number of domains is not 2 or 3")

        return tabulate(rows, headers=headers, tablefmt="github")

    def __getitem__(self, key: tuple[Hashable | None, ...]) -> float:
        """Return the score for a symbol tuple, with fallback for missing keys.

        Args:
            key: Symbol tuple to look up.

        Returns:
            The score, or the minimum matching score if key is not found.
        """
        if key not in self.scores:
            potential = []
            for domain, state in enumerate(key):
                for entry, score in self.scores.items():
                    if entry[domain] == state:
                        potential.append(score)

            if not potential:
                return min(self.scores.values())
            return min(potential)

        return self.scores[key]
