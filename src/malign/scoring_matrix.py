"""
Module for scoring matrices.
"""

# Import Python standard libraries
from typing import Dict, Hashable, List, Optional, Tuple
import copy
import itertools
import json

# Import 3rd-party libraries; `enable_iterative_imputer` is necessary for the one in `sklearn.impute`
from sklearn.experimental import (
    enable_iterative_imputer,
)  # pylint: disable=unused-import
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.impute import IterativeImputer, SimpleImputer
from sklearn.linear_model import BayesianRidge
from sklearn.neighbors import KNeighborsRegressor
from sklearn.tree import DecisionTreeRegressor
from tabulate import tabulate
import numpy as np


# TODO: study different method for loading


class ScoringMatrix:
    """
    Class for sequence alignment scoring matrices.

    A scoring matrix is implemented to work as a Python dictionary and, in
    fact, the alignment methods should be written as expecting any
    data structure that works like a dictionary, via the `__get__()`
    method.

    The class, however, provides a number of facilities for dealing with the
    multidimensionality of the matrices for multiple alignment on
    multiple domains, as intended by the library. Those are mainly:

      - allowing and mandating a unique domain for each sequence, even in
        cases where the domains are equal
      - allow asymmetric scoring value
      - treat gaps as full, normal "symbols"

    The class especially facilitates dealing with sparse matrices/vectors and
    with submatrices, both for querying and for building.

    The `gap` symbol is a mandatory property and must be shared by all
    domains. By design, exclusive gap alignment vectors (that is, keys
    composed only of gaps in all the provided domains) will have a value of
    0.0.

    For incomplete vectors relating to missing or non-informed domains, these
    must be specified by using a `None` value in the key.
    """

    def __init__(
        self,
        scores: Dict[Tuple[Hashable, ...], float] = None,
        domains: Optional[List[List[Hashable]]] = None,
        gap: Hashable = "-",
        impute_method: Optional[str] = "mean",
    ):
        """
        Initialize a scoring matrix.

        @param scores: An optional scoring dictionary, with tuples of hashable elements, usually
            strings, as keys (respecting the order that will be used when using the matrix) and
            floating points as values. If not provided, a blank ScoringMatrix will be initialized,
            and all other initialization arguments will be disregarded.
        @param domains: A list of lists of hashable elements, with the set of symbols ("domain") that
            composes each domain that will be used for alignment. If not provided, it will be
            computed from the `scores`, assuming that all symbols are used at least once; the
            only symbol automatically added to each domain, if not found in `scores`, is the `gap`.
            Note that by definition, as it is treated as a normal value, the gap symbol must be
            part of each domain unless it is explicitly assumed that the domain does not allow gaps.
        @param gap: The hashable element to use as a gap symbol. Defaults to `"-"`.
        @param impute_method: An optional string with the name of the method used to fill empty
            elements in the matrix, using values imputed from the ones that are provided. If not
            set (i.e., if equal to `None`), the matrix will not be filled. Choices are `mean`,
            `median`, `decision_tree`, `extra_trees`, `k_neighbors`, and `bayesian_ridge`.
            Defaults to `mean`.
        """

        # Set properties provided by the user, so that the matrix can be used even empty
        self.impute_method = impute_method
        self.gap = gap

        # Instantiate properties that will be initialized later
        # TODO: replace `._dr` with a method reading domains
        self.num_domains: int = -1
        self.domains: Optional[List[List[Hashable]]] = domains
        self.score = scores
        self._dr: Optional[Tuple[int, ...]] = None

        # If `scores` is not provided, we initialize an empty scorer (which will likely load from disk).
        if scores:
            # Extract the domains (or build them, if they were not provided) and
            # make sure they are sorted and comparable.
            self._init_domains(domains, scores)

            # Initialize from the scores
            self._init_from_scores(scores)

    def _init_domains(
        self,
        domains: Optional[List[List[Hashable]]],
        scores: Dict[Tuple[Hashable, ...], float],
    ):
        """
        Internal method for initializing domains from user-provided arguments.

        @param domains: The `domains` argument passed to `.__init__()`.
        @param scores: The `scores` argument passed to `.__init__()`.
        """

        if domains is not None:
            self.domains = domains
        else:
            # Collect domains from scores. Note that the gap symbol is the
            # only one automatically added to each domain, and that
            # `None`s are discarded as the user or the library might be
            # passing sub-matrices as well (as when creating an
            # identity matrix, when it makes more sense to provide the
            # sub-matrix identities.
            self.domains = [
                set([symbol for symbol in domain if symbol] + [self.gap])
                for domain in zip(*scores.keys())
            ]

        # Organize the domains; as this is applied to both user-provide
        # and code-derived lists, we apply `set` and `sorted` in all
        # cases to simplify the code
        self.domains = [sorted(set(domain)) for domain in self.domains]

    def _init_from_scores(self, scores: Dict[Tuple[Hashable, ...], float]):
        """
        Internal method for initializing from user-provided scores.

        The method assumes that domains have already been initialized by an internal
        call to `._init_domains()`.

        @param scores: The `scores` argument passed to `.__init__()`.
        """

        # Make sure all `scores` report the same number of domains, store
        # the number of domains and an internal variable with the range (which
        # is used multiple times; the one facing the user is `self.num_domains`)
        # TODO: remove "domain range" (._dr) or, at least, expand the name
        self.num_domains = len(self.domains)
        self._dr = tuple(range(self.num_domains))

        key_lens = {len(key) for key in scores}
        if len(key_lens) > 1:
            raise ValueError("Different domain-lengths in `scores`.")

        # Make sure the domains contain all symbols used in `scores`; while
        # this would not be strictly necessary when the domains are collected
        # automatically, it is still good to keep a simpler loop for that
        # in all circumstances
        scores_domains = list(zip(*scores.keys()))
        found = [
            all([symbol in ref_domain for symbol in scores_domain if symbol])
            for scores_domain, ref_domain in zip(scores_domains, self.domains)
        ]
        if not all(found):
            raise ValueError("`scores` has symbols not in domains.")

        # Store a copy of the `scores` and set (or override) the value of the
        # all-gap vector
        self.scores = scores.copy()
        self.scores[tuple([self.gap] * self.num_domains)] = 0.0

        # Fill the matrix with the appropriate method if requested.
        if self.impute_method:
            if self.impute_method not in [
                "mean",
                "median",
                "decision_tree",
                "extra_trees",
                "k_neighbors",
                "bayesian_ridge",
            ]:
                raise ValueError(f"Unknown imputation method: {self.impute_method}.")

            self._fill_matrix()

    def _fill_matrix(self):
        """
        Internal function for filling a matrix with an imputation method.
        """

        # We perform matrix imputation on multi-hot vectors, with one true
        # position for each domain. `None`, as used for submatrices, is always
        # mapped to the first site of its domain; the gap symbol will be part
        # of the sorted list of symbols, if found. In order to do that, we
        # build an initial list of all domain/values, so we can fill the
        # vectors. We call it `encoder` due to analogy with data analysis
        # pipelines, even though it is not strictly an encoder but an auxiliary
        # data structure.
        encoder = list(
            itertools.chain.from_iterable(
                [
                    [(domain_idx, value) for value in [None] + domain_values]
                    for domain_idx, domain_values in enumerate(self.domains)
                ]
            )
        )

        # Fill the matrix for imputation for each possible combination
        domains_with_none = [[None] + domain for domain in self.domains]
        train_matrix = []
        imp_matrix = []
        for cat_vector in itertools.product(*domains_with_none):
            # We need to add `None`s in order to account for submatrices,
            # but also need to make sure that at least two non-`None`s are
            # found (we need at least to sequences to have an alignment)
            if len([val for val in cat_vector if val]) < 2:
                continue

            # Build the multi-hot vector
            mh_vector = [cat_vector[element[0]] == element[1] for element in encoder]
            score = self.scores.get(cat_vector)
            if score is not None:
                train_matrix.append(mh_vector + [score])
            else:
                imp_matrix.append(mh_vector + [np.nan])

        # If the imputation matrix is empty, there is nothing to do
        if not imp_matrix:
            return

        # Instantiate the imputer, train it, and fill the matrix
        # TODO: allow setting more options when initializing, such num estimators and random state
        if self.impute_method == "decision_tree":
            estimator = DecisionTreeRegressor(max_features="sqrt", random_state=0)
            imputer = IterativeImputer(random_state=0, estimator=estimator)
        elif self.impute_method == "extra_trees":
            estimator = ExtraTreesRegressor(n_estimators=10, random_state=0)
            imputer = IterativeImputer(random_state=0, estimator=estimator)
        elif self.impute_method == "k_neighbors":
            estimator = KNeighborsRegressor(n_neighbors=15)
            imputer = IterativeImputer(random_state=0, estimator=estimator)
        elif self.impute_method == "bayesian_ridge":
            estimator = BayesianRidge()
            imputer = IterativeImputer(random_state=0, estimator=estimator)
        elif self.impute_method in ["mean", "median"]:
            imputer = SimpleImputer(missing_values=np.nan, strategy=self.impute_method)
        else:
            # Currently default to `mean`
            imputer = SimpleImputer(missing_values=np.nan, strategy="mean")

        imputer.fit(train_matrix)
        trans_matrix = imputer.transform(imp_matrix)

        # Reverse from the transformed matrix and store the value
        for row in trans_matrix:
            mh_vector, value = row[:-1], row[-1]

            # `encoder[idx]` is a tuple (position, character), but we only need
            # the character as we can trust the code not to give us more than
            # one character per position
            cat_vector = [
                encoder[idx][1] for idx, value in enumerate(mh_vector) if value
            ]

            self.scores[tuple(cat_vector)] = value

    # TODO: allow to return a ScoringMatrix, so we have `ScoringMatrix().load(filename)`
    # TODO: make sure it works with types other than strings, at least integers
    # TODO: better solution than JSON for serializing? perhaps just pickle?
    def load(self, filename: str) -> None:
        """
        Internal function for loading a serialized matrix from file.

        @param filename: Path to the file holding the matrix to be loaded.
        """

        with open(filename, encoding="utf-8") as json_handler:
            serial_data = json.load(json_handler)

            self.gap = serial_data["gap"]
            self._dr = tuple(serial_data["domain_range"])
            self.domains = serial_data["domains"].copy()

            # Make sure values are floats
            # TODO: should "NULL" be "0000000000"? "NULL" could be a symbol...
            self.scores = {
                tuple(
                    [
                        None if sub_key == "NULL" else sub_key
                        for sub_key in key.split(" / ")
                    ]
                ): float(value)
                for key, value in serial_data["scores"].items()
            }

    def save(self, filename: str) -> None:
        """
        Serialize a matrix to disk.

        Note that, due limits in serializing with JSON, which do not allow
        lists as keys, the `" / "` symbol used for serialization (note the
        spaces) is not allowed.

        @param filename: Path to the file where to write the serialized matrix.
        """

        # Make sure the reserved symbol is not used
        if any([" / " in domain for domain in self.domains]):
            raise ValueError("At least one domain uses reserved symbols.")

        # Make a copy of `self.scores` replacing `None`s, which cannot be part
        # of the key as per the JSON standard, with a custom element, which
        # must be mapped back when loading.
        # As we only operate on strings, to facilitate human edition/inspection
        # we map the `None`s to numeric zero value. Note that we need our
        # custom sorting function as we cannot otherwise sort keys with Nones
        # and strings.
        def _allow_none_key(obj):
            return tuple(["0000000000" if k is None else k for k in obj])

        _scores = {
            " / ".join(["NULL" if k is None else k for k in key]): self.scores[key]
            for key in sorted(self.scores, key=_allow_none_key)
        }

        # Open handler, build serialized data, and write to disk
        with open(filename, "w", encoding="utf-8") as json_handler:
            serial_data = {
                "domains": self.domains,
                "gap": self.gap,
                "domain_range": list(self._dr),
                "scores": _scores,
            }

            json.dump(serial_data, json_handler, indent=4, ensure_ascii=False)

    # TODO: should cache? should be internal?
    # TODO: cannot annotate that it returns a dict with ScoringMatrix as values
    def compute_submatrices(self, domains: List[Tuple[int, ...]]) -> Dict:
        """
        Compute submatrices from a collection of domains.

        Sub-matrices are useful, when compared to the possibility of addressing
        a normal matrix with `None` values in keys, in simplifying computation
        and debugging, as well as in some methods for smoothing.

        @param domains: A list of domains for which to compute submatrices, as tuples.
        @return: A dictionary with domain-tuples as keys and ScoringMatrices as values.
        """

        sub_matrix_scores: Dict[Tuple[Hashable, ...], ScoringMatrix] = {}

        for sub_domain in domains:
            # Trigger the computation/completion of the domain
            # TODO: extend documentation to note that this assumes the matrix has been filled
            # TODO: maybe some flag? check if filled indeed?

            # Collect all scores for the given domain
            sub_scores: Dict[Tuple[Hashable, ...], float] = {}
            for key, value in self.scores.items():
                # Make sure that all entries of domain exist...
                check = [key[idx] is not None for idx in sub_domain]
                # ...and all others don't
                check += [key[idx] is None for idx in self._dr if idx not in sub_domain]

                if all(check):
                    sub_key = tuple([k for k in key if k is not None])
                    sub_scores[sub_key] = value

            # Create the new ScoringMatrix
            sub_matrix_scores[sub_domain] = ScoringMatrix(sub_scores)

        return sub_matrix_scores

    # TODO: use __copy__ and __deepcopy__ properly
    # TODO: cannot annotate it returns a scoring matrix
    def copy(self):
        """
        Return a copy of the current matrix object.
        """

        # TODO: can we perform the copy manually?
        return copy.deepcopy(self)

    # TODO: currently only working for 2 or 3 domains
    # TODO: use more options from tabulate
    def tabulate(self) -> str:
        """
        Build a string with a tabulated representation of the matrix.

        @return: A string with a tabular representation of the matrix, intended for human
            inspection.
        """

        rows = []
        if self.num_domains == 2:
            for symbol_a in self.domains[0]:
                row = [symbol_a] + [
                    self.scores[symbol_a, symbol_b] for symbol_b in self.domains[1]
                ]
                rows.append(row)

            headers = [""] + list(self.domains[1])

        elif self.num_domains == 3:
            rows = []
            for symbol_a in self.domains[0]:
                row = [symbol_a] + [
                    self.scores.get((symbol_a, symbol_b, symbol_c), "-")
                    for symbol_b, symbol_c in itertools.product(
                        self.domains[1], self.domains[2]
                    )
                ]
                rows.append(row)

            headers = [""] + [
                "/" + "/".join(sub_key)
                for sub_key in itertools.product(self.domains[1], self.domains[2])
            ]
        else:
            raise ValueError("number of domains is not 2 or 3")

        return tabulate(rows, headers=headers, tablefmt="github")

    # TODO: decide/inform on what to do when a sub-domain is requested and it has
    #       not been computed.
    def __getitem__(self, key: Tuple[Optional[Hashable], ...]) -> float:
        """
        Return the score associated with a tuple of alignments per domain.

        Sub-domains can be specified by passing a `None` to the positions the key
        does not apply, e.g. `("a", "b", None)`.

        @param key: The element in the multidimensional matrix to be queried.
        @return: The value associated with the requested element.
        """

        return self.scores[key]

    def __setitem__(self, key: Tuple[Hashable, ...], value: float) -> None:
        """
        Set the value of a key in the multidimensional matrix.

        @param key: The element in the multidimensional matrix whose value will be set or changed.
        @param value: The value to which set the element.
        """

        matches = [k in domain for k, domain in zip(key, self.domains) if k is not None]
        if not all(matches):
            raise ValueError("`key` uses symbol(s) not in domain.")

        self.scores[key] = value
