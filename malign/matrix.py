"""
Module for scoring matrices.
"""

# Import Python standard libraries
from collections import defaultdict
import itertools
import json
import copy

# Import 3rd-party libraries
import numpy as np
from tabulate import tabulate

from . import utils

# TODO: add methods for loading and storing the matrices
# TODO: add auxiliary function to build a scoring matrix in a single pass from subs
# TODO: add new method for inference using something like impurity/entropy
# TODO: implement a __get__ method, related to __call__

# TODO: `scores` can be empty if submatrices are provided
# TODO: what if submatrices are not complete?
class ScoringMatrix:
    """
    Class for sequence alignment scoring matrices.

    A scoring matrix is implemented to work as a Python dictionary and, in
    fact, it should be allowed to use a plain dictionary instead of a
    ScoringMatrix. The class provides facilities for dealing with the
    multidimensionality of the matrices for multiple alignmentes used in
    this library, that is:

      - a unique alphabet for each sequence
      - asymmetric information
      - gaps are full, normal "symbols"
    """

    def __init__(self, scores=None, filename=None, **kwargs):
        """
        Initialize a scoring matrix.

        Parameters
        ==========

        scores : dict
            A scoring dictionary, with tuples of strings as keys (respecting
            the order that will be used when using the matrix) and floating
            points as values. Exclusive-gap alignments (that is, keys
            composed only of gaps) will be overriden to a value of 0.0 if
            provided. Either `scores` or `filename` must be provided when initializating
            a matrix.
        filename : str
            A string with the path to a serialized matrix in JSON format. Either
            `scores` or `filename` must be provided when initializing a matrix.
        sub_matrices : dict
            A dictionary with domains (tuples of ints) as scorers and
            ScoringMatrices as values. Assumes the sub-matrices are
            complete.
        alphabets : list of lists of strings
            A list of lists of strings, with the alphabets of each domain
            that will be used for alignment. If provided, the alphabet lists
            will be sorted. If not provided, it will be
            computed from the `scores`, assuming that all symbols are used
            at least once. Note that, by definition, the gap symbol must be
            part of the alphabet (unless it is explicitly assumed that the
            domain does not allow gaps).
        gap : str
            The string to use as a gap symbol. Defaults to `-`.
        fill_method : string or None
            The method used to fill empty elements in the matrix, using
            values computed from the ones provided. If set to `None`, the
            matrix will not be filled. Choices are `"standard"` and `"distance"`.
            Defaults to the `"standard"` method.
        """

        # TODO: check precedence/interaction filename/scorer/submatirces

        # Fill the matrix with the appropriate method if requested
        self._fill_method = kwargs.get("fill_method", "standard")

        # If a filename was provided, load a serialized matrix; otherwise, initialize
        # from user provided `scores`
        if filename:
            self._load(filename)
        else:
            # Store values
            self.gap = kwargs.get("gap", "-")

            # Collect submatrices, if they were provided
            sub_matrices = kwargs.get("sub_matrices", {})

            # Extract the alphabets, if they were not provided; during extraction,
            # `None`s are removed as the user or the library might be passing
            # sub-matrices as well (as when creating an identity matrix, when it
            # makes more sense to provide the sub-matrix identities)
            # TODO: check if the alphabets agree with the provided scores
            # TODO: check if the alphabets are list of lists, if provided
            self.alphabets = kwargs.get("alphabets", None)

            if self.alphabets:
                self.alphabets = [
                    tuple(sorted(alphabet)) for alphabet in self.alphabets
                ]
            else:
                # Collect alphabet from scores; if no `scores` were provided but
                # only `sub_matrices`, we need to initialize `self.alphabets` to the
                # length we can derive from the `sub_matrices` keys
                if scores:
                    self.alphabets = [alphabet for alphabet in zip(*scores.keys())]
                else:
                    num_alphabets = max([max(domain) for domain in sub_matrices])
                    self.alphabets = [[] for _ in range(num_alphabets + 1)]

                # Extend with alphabets from submatrices, if they were provided
                for sub_matrix, obj in sub_matrices.items():
                    for dmn, alphabet in zip(sub_matrix, obj.alphabets):
                        self.alphabets[dmn] += alphabet

                # Make sorted lists, making sure the gap is there
                self.alphabets = [
                    tuple(
                        sorted(
                            set(
                                [symbol for symbol in alphabet if symbol is not None]
                                + [self.gap]
                            )
                        )
                    )
                    for alphabet in self.alphabets
                ]

            # Initialize from the scores
            self._init_from_scores(scores, sub_matrices)

    def _init_from_scores(self, scores, sub_matrices):
        """
        Internal function for initializing from user-provided scores.
        """

        # Make sure all `scores` report the same number of domains, store
        # the number of domains, an set (or override, if necessary) the
        # value of exclusive-gap sites.
        # TODO: check if in line with submatrices
        if scores:
            domains = set([len(key) for key in scores])
            if len(domains) > 1:
                raise ValueError("Different domain-lengths in `scores`.")

            # Copy `domains` and define the global, full `domain range`
            # NOTE: the comma after `self.domains` is used to unpack the
            # `domains` set, which at this point we know contains a single
            # item
            self.domains, = domains
        else:
            self.domains = 1 + max([max(domain) for domain in sub_matrices])

        #
        self._dr = tuple(range(self.domains))
        self.scores = scores.copy()
        self.scores[tuple([self.gap] * self.domains)] = 0.0

        # Add submatrices, if provided
        for sub_domain, matrix in sub_matrices.items():
            # We need to index the submatrix with its `._dr`, as the indexes
            # we need are the ones in this matrix being initialized. Note
            # how we keep a single dictionary, filling missing spots with
            # `None`
            for sub_key, score in matrix.scores.items():
                mapper = {sub_idx: idx for sub_idx, idx in zip(sub_domain, matrix._dr)}
                new_key = [
                    sub_key[idx] if idx is not None else None
                    for idx in [mapper.get(idx, None) for idx in self._dr]
                ]

                self.scores[tuple(new_key)] = score

        # Fill the matrix with the appropriate method if requested
        if self._fill_method:
            if self._fill_method in ["standard", "distance"]:
                self._fill_full_matrix(self._fill_method)
            else:
                raise ValueError("Unknown filling method.")

    # TODO: currently disregarding the `method`, as there is a single one
    # TODO: describe how the standard method can be considered a kind of MLE
    # TODO: change method to first look for the closest match, especially when
    #       we have submatrices
    # TODO: check A Fill Estimation Algorithm for Sparse Matrices and Tensors in Blocked Formats
    # TODO: another method could be a weighted avarage from the distance, which
    #       can be computed in terms of neighbouts, involving the entire matrix
    def _fill_full_matrix(self, method):
        """
        Internal function for filling a matrix if there are missing values.

        Parameters
        ==========

        method : string
            The method to be used for filling the matrix. Choices are
            `"standard"`.
        """

        # TODO: describe distance method
        if self._fill_method == "distance":
            score_cache = {}
            for cur_key in itertools.product(*self.alphabets):
                if cur_key not in self.scores:

                    # Collect all weighted distances
                    x = []
                    for ref_key, score in self.scores.items():
                        match = sum(
                            [
                                symbol_key == symbol_ref
                                for symbol_key, symbol_ref in zip(ref_key, cur_key)
                            ]
                        )
                        x += [score] * match

                    # TODO: allow mean
                    score_cache[cur_key] = np.mean(x)

            # Update an return
            # TODO: check if it fails when we have no info (extreme case), just
            # setting a mean value should be enough (this could be for all methods)
            self.scores.update(score_cache)
            return

        # Compute the expected size of the scorer, so we know if there are
        # keys missing indeed; we check for `all` to exclude submatrices
        # that have `None` values
        # TODO: deal with the new Nones (should just filter?), also below
        expected_size = np.prod([len(alphabet) for alphabet in self.alphabets])
        if len([key for key in self.scores if all(key)]) == expected_size:
            return

        # For each domain, we first collect all keys without considering
        # the symbol in the domain itself, and fill any missing spot with
        # the mean value. No difference is made in terms of gaps.
        # Example: if we have ('a', 'A', '1')=10 and ('b', 'A', '1')=0, but
        # no ('c', 'A', '1'), this will set the latter to the mean value
        # of 5, as coming from (?, 'A', '1')
        sub_scores = {}
        for domain_idx in self._dr:
            # Collect all sub-keys and values
            sub_scores[domain_idx] = defaultdict(list)
            for key, score in self.scores.items():
                # Skip full gap
                if all([v == self.gap for v in key]):
                    continue

                sub_key = tuple(
                    [value for idx, value in enumerate(key) if idx != domain_idx]
                )

                # We might have `None`s due to submatrices
                if not None in sub_key:
                    sub_scores[domain_idx][sub_key].append(score)

        # Set the new score when possible, from the mean of the
        # available `sub_scores`. We cache the new values in order to only
        # apply them after the loop, not changing `self.scores` in-place.
        score_cache = {}
        for key in itertools.product(*self.alphabets):
            if key not in self.scores:
                # Collect all sub-scores
                all_sub_scores = []
                for domain_idx in self._dr:
                    sub_key = key[:domain_idx] + key[domain_idx + 1 :]
                    all_sub_scores += sub_scores[domain_idx].get(sub_key, [])

                # Cache the new value if possible
                if any(all_sub_scores):
                    # TODO: allow other methods, such as mean/median
                    score_cache[key] = np.percentile(all_sub_scores, 50)

        # Update with the new values
        self.scores.update(score_cache)

        # If there are still keys missing, collect the mean sub_score
        # by symbol in each domain and set to this value -- this is similar
        # to what we do above, but much less precise as we take all symbols
        # independently
        if len([key for key in self.scores if all(key)]) == expected_size:
            return

        symbol_score = defaultdict(lambda: defaultdict(list))
        for key, score in self.scores.items():
            for domain_idx in self._dr:
                symbol_score[domain_idx][key[domain_idx]].append(score)

        for domain_idx in self._dr:
            symbol_score[domain_idx] = {
                symbol: np.mean(scores)
                for symbol, scores in symbol_score[domain_idx].items()
            }

        for key in itertools.product(*self.alphabets):
            if key not in self.scores:
                # If the alphabet passed by the user has symbols not in the scorer,
                # we will have a KeyError in symbol_score[domain_idx][symbol]; the
                # `.get()` will default towards the `gap`, which is always necessary
                # TODO: find a better solution
                self.scores[key] = np.mean(
                    [
                        symbol_score[domain_idx].get(
                            symbol, symbol_score[domain_idx][self.gap]
                        )
                        for domain_idx, symbol in enumerate(key)
                    ]
                )

    # TODO: if the alphabets are the same, we can save a lot of computation,
    #       by caching, such as in identity matrices
    def _fill_domain(self, domain):
        """
        Internal function for filling a sub-domain.

        The function assumes the full domain matrix has already been
        filled.

        Parameters
        ==========

        domain : tuple of ints
            Tuble of the sub-domain to be filled.
        """

        # Obtain the alphabets of the sub-domain
        sub_alphabets = [
            self.alphabets[d_idx] if d_idx in domain else [None] for d_idx in self._dr
        ]

        # Fill keys
        # TODO: investigate better subsetting, loops to unroll
        # TODO: also what should not be collected like full gaps
        for sub_key in itertools.product(*sub_alphabets):
            if sub_key not in self.scores:
                sub_scores = []
                for key, score in self.scores.items():
                    if all(
                        [
                            key[idx] == v
                            for idx, v in enumerate(sub_key)
                            if v is not None
                        ]
                    ):
                        sub_scores.append(score)

                # Add adjusted value
                # TODO: allow other manipulations? percentile perhaps?
                self.scores[sub_key] = max(sub_scores)

    def _load(self, filename):
        with open(filename) as json_handler:
            serial_data = json.load(json_handler)

            self.gap = serial_data["gap"]
            self._dr = tuple(serial_data["domain_range"])
            self.alphabets = [tuple(alphabet) for alphabet in serial_data["alphabets"]]

            self.scores = {
                tuple([None if k == "NULL" else k for k in key.split(" / ")]): value
                for key, value in serial_data["scores"].items()
            }

    # TODO: note issues about " / " if in alphabets
    def save(self, filename):
        # Make a copy of `self.scores` replacing `None`s, which cannot be part
        # of the key as per the JSON standard, with a custom element, which
        # must be mapped back when loading.
        # As we only operate on strings, to facilitate human edition/inspection
        # we map the `None`s to numeric zero value. Note that we need our custom
        # sorting function as we cannot otherwise sort keys with Nones and strings.
        def _allow_none_key(obj):
            return tuple(["0000000000" if k is None else k for k in obj])

        _scores = {
            " / ".join(["NULL" if k is None else k for k in key]): self.scores[key]
            for key in sorted(self.scores, key=_allow_none_key)
        }

        # Build serialized data and write to disk
        serial_data = {
            "alphabets": self.alphabets,
            "gap": self.gap,
            "domain_range": list(self._dr),
            "scores": _scores,
        }
        with open(filename, "w") as json_handler:
            json.dump(serial_data, json_handler, indent=4, ensure_ascii=False)

    # compute the submatrices as a dictionary of ScoringMatrix
    def compute_submatrices(self, domains):
        sub_matrix = {}

        for domain in domains:
            # Trigger the computation/completion of the domain
            self._fill_domain(domain)

            # Collect all scores for the given domain
            sub_scores = {}
            for key, value in self.scores.items():
                # Make sure all entries of domain exist...
                check = [key[idx] is not None for idx in domain]
                # ...and all others don't
                check += [key[idx] is None for idx in self._dr if idx not in domain]

                if all(check):
                    sub_key = tuple([k for k in key if k is not None])
                    sub_scores[sub_key] = value

            # Create the ScoringMatrix
            sub_matrix[domain] = ScoringMatrix(sub_scores)

        return sub_matrix

    def copy(self):
        return copy.deepcopy(self)

    # TODO: currently only working for 2 or 3 domains
    # TODO: use more options from tabulate
    # TODO: check about identity matrix? from alphabets with a method?
    def tabulate(self):
        rows = []
        if len(self._dr) == 2:

            for symbol_a in self.alphabets[0]:
                row = [symbol_a] + [
                    self.scores[symbol_a, symbol_b] for symbol_b in self.alphabets[1]
                ]
                rows.append(row)

            headers = [""] + list(self.alphabets[1])

        elif len(self._dr) == 3:
            rows = []
            for symbol_a in self.alphabets[0]:
                row = [symbol_a] + [
                    self.scores.get((symbol_a, symbol_b, symbol_c), "-")
                    for symbol_b, symbol_c in itertools.product(
                        self.alphabets[1], self.alphabets[2]
                    )
                ]
                rows.append(row)

            headers = [""] + [
                "/" + "/".join(sub_key)
                for sub_key in itertools.product(self.alphabets[1], self.alphabets[2])
            ]
        else:
            raise ValueError("number of domains is not two or 3")

        return tabulate(rows, headers=headers, tablefmt="github")

    def __getitem__(self, key):
        """
        Return the score associated with a tuple of alignments per domain.

        Sub-domains can be specified by passing a `None` to the positions the key
        does not apply.

        If a sub-domain is requested and it has not been computed so far,
        it will be inferred and cached for reuse.

        Parameters
        ==========

        key : tuple of strings
            The element in the multidimensional matrix to be queried.

        Returns
        =======

        value : float
            The value associated with the element.
        """

        # If there are `None`s in the key, it is a subdomain, be sure to check
        # if it needs to be computed

        # If there are `None`s in the key (as observed from the domain), and the key
        # is missing, it means the submatrix was not provided and has not been computed
        # yet; we do it now, caching the results.
        if None in key:
            domain = tuple([idx for idx, value in enumerate(key) if value is not None])

            if key not in self.scores:
                self._fill_domain(domain)

        return self.scores[tuple(key)]

    # TODO: run checks, alphabet, etc.
    def __setitem__(self, key, value):
        self.scores[tuple(key)] = value
