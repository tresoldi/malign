"""
Module for scoring matrices.
"""

# Import Python standard libraries
from collections import defaultdict
import copy
import itertools
import json

# Import 3rd-party libraries
import numpy as np
from tabulate import tabulate


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

    The object especially facilitates dealing with sparse matrices/vectors and
    with submatrices (i.e., not full domains), both for querying and for
    building.
    """

    def __init__(self, scores=None, **kwargs):
        """
        Initialize a scoring matrix.

        Either the `filename` or a combination of `scores` and `submatrices` must be
        provided. At least one of the values of `scores` and `submatrices` must be
        provided.

        Parameters
        ==========

        scores : dict
            A scoring dictionary, with tuples of strings as keys (respecting
            the order that will be used when using the matrix) and floating
            points as values. Exclusive-gap alignments (that is, keys
            composed only of gaps) will be overriden to a value of 0.0 if
            provided. Subdomains are indicated by `None` valus in the keys.
            Defaults to `None`.
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

        # Extract `scores`, `filename`, and `submatrices`, if provided
        filename = kwargs.get("filename", None)
        sub_matrices = kwargs.get("sub_matrices", {})

        # If a filename was provided, load a serialized matrix; otherwise, initialize
        # from user provided `scores`
        if filename:
            self._load(filename)
        else:
            # Extract additional values; note that these, as everything else, are
            # not considered when loading from disk

            # Store values
            self._fill_method = kwargs.get("fill_method", "standard")
            self.gap = kwargs.get("gap", "-")

            # Extract the alphabets or build them, if they were not provided;
            # during extraction, `None`s are removed as the user or the library might
            # be passing sub-matrices as well (as when creating an identity matrix,
            # when it makes more sense to provide the sub-matrix identities)
            self.alphabets = kwargs.get("alphabets", None)

            if self.alphabets:
                # Organize provided alphabets
                self.alphabets = [sorted(alphabet) for alphabet in self.alphabets]

                # Make sure the alphabets contain all symbols used in `scores`
                # TODO: add check for submatrices' symbols, if provided
                scores_alphabets = list(zip(*scores.keys()))
                available = [
                    all([symbol in ref_alphabet for symbol in scores_alphabet])
                    for scores_alphabet, ref_alphabet in zip(
                        scores_alphabets, self.alphabets
                    )
                ]
                if not all(available):
                    raise ValueError("`scores` has symbols not in alphabets.")

            else:
                # Collect alphabet from scores; if no `scores` were provided but
                # only `sub_matrices`, we need to initialize `self.alphabets` to the
                # length we can derive from the `sub_matrices` keys
                if scores:
                    self.alphabets = [
                        list(alphabet) for alphabet in zip(*scores.keys())
                    ]
                else:
                    num_alphabets = max([max(domain) for domain in sub_matrices])
                    self.alphabets = [[] for _ in range(num_alphabets + 1)]

                # Extend with alphabets from submatrices, if they were provided
                for sub_matrix, obj in sub_matrices.items():
                    for dmn, alphabet in zip(sub_matrix, obj.alphabets):
                        self.alphabets[dmn] += alphabet

                # Make sorted lists, making sure the gap is there
                self.alphabets = [
                    sorted(
                        set(
                            [symbol for symbol in alphabet if symbol is not None]
                            + [self.gap]
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
        # the number of domains
        if scores:
            domains = {len(key) for key in scores}
            if len(domains) > 1:
                raise ValueError("Different domain-lengths in `scores`.")

            # Copy `domains` and define the global, full `domain range`
            # NOTE: the comma after `self.domains` is used to unpack the
            # `domains` set, which at this point we know contains a single
            # item
            self.domains, = domains
        else:
            self.domains = 1 + max([max(domain) for domain in sub_matrices])

        # Store a copy of the `scores` and cache the `domain_range` as an internal
        # variable (as it will be used repeated times -- the one facing the user
        # is `self.domains`)
        self.scores = scores.copy()
        self._dr = tuple(range(self.domains))

        # Set (or override) the value of all-gap vector
        self.scores[tuple([self.gap] * self.domains)] = 0.0

        # Add submatrices, if provided
        for sub_domain, matrix in sub_matrices.items():
            # We need to index the submatrix with its `._dr`, as the indexes
            # we need are the ones in this matrix being initialized. Note
            # how we keep a single dictionary, filling missing spots with `None`
            # pylint: disable=protected-access
            mapper = dict(zip(sub_domain, matrix._dr))
            for sub_key, score in matrix.scores.items():
                # sub_key ('b', 'X') -> new_key ['b', 'X', None]
                # sub_key ('c', '-') -> new_key ['c', '-', None]
                # sub_key ('a', '-') -> new_key ['a', None, '-']
                new_key = [
                    sub_key[idx] if idx is not None else None
                    for idx in [mapper.get(idx, None) for idx in self._dr]
                ]

                # Make sure all values are floats
                self.scores[tuple(new_key)] = float(score)

        # Fill the matrix with the appropriate method if requested. Note that it is
        # not easy to verify beforehand if the matrix needs to be filled, as the
        # combination of submatrices and different domains implies that the
        # expected size is not necessarily the product of the alphabets. For making
        # the code simpler to follow, no check is performed beforehand.
        if self._fill_method:
            if self._fill_method == "standard":
                self._fill_matrix_standard()
            elif self._fill_method == "distance":
                self._fill_matrix_distance()
            else:
                raise ValueError("Unknown filling method.")

            # Run the fallback for cells we could not fit, if any
            self._fill_matrix_fallback()

    def _fill_matrix_standard(self):
        """
        Internal function for filling a matrix with the `standard` method.
        """

        # For each domain, we first collect all keys without considering
        # the symbol in the domain itself, and fill any missing spot with
        # the some descriptive value (mean, median, percentile, etc.).
        # No difference is made in terms of treatment of gaps.
        # Example: if we have ('a', 'A', '1')==10 and ('b', 'A', '1')==0, but
        # no ('c', 'A', '1'), this will set the latter to the mean value
        # of 5, as coming from (?, 'A', '1')
        sub_scores = {}
        for domain_idx in self._dr:
            # Collect all sub-keys and values
            sub_scores[domain_idx] = defaultdict(list)
            for key, score in self.scores.items():
                # Don't account for the full gap vector
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
                    # TODO: allow other methods, such as percentile/median
                    score_cache[key] = np.mean(all_sub_scores)

        # Update with the new values
        self.scores.update(score_cache)

    def _fill_matrix_distance(self):
        """
        Internal function for filling a matrix with the `distance` method.
        """

        score_cache = {}
        for cur_key in itertools.product(*self.alphabets):
            if cur_key not in self.scores:

                _scores = []
                for ref_key, score in self.scores.items():
                    match_sites = sum(
                        [
                            symbol_key == symbol_ref
                            for symbol_key, symbol_ref in zip(ref_key, cur_key)
                        ]
                    )
                    _scores += [score] * match_sites

                # TODO: allow measures other than mean
                score_cache[cur_key] = np.mean(_scores)

        # Update with the new values
        self.scores.update(score_cache)

    def _fill_matrix_fallback(self):
        """
        Internal function for last resource on matrix filling.

        This function will be used when there are still empty cells in cell and
        the method for filling it was unable to account for. As a last resource,
        it is intended as a common method.
        """

        # Holder for temporary scores
        symbol_score = defaultdict(list)
        gap_scores = []

        # Collect a dictionary of all scores for all symbols in each domain, plus
        # the mean score overall for gaps
        # NOTE: this is collecting also full-gap alignment
        for key, score in self.scores.items():
            for domain_idx in self._dr:
                # Collect cell score
                symbol_score[domain_idx, key[domain_idx]].append(score)

                # Collect score if a gap
                if key[domain_idx] == self.gap:
                    gap_scores.append(score)

        symbol_score = {key: np.mean(scores) for key, scores in symbol_score.items()}
        gap_score = np.mean(gap_scores)

        for key in itertools.product(*self.alphabets):
            if key not in self.scores:
                # If the alphabet passed by the user has symbols not in the scorer,
                # we will have a KeyError in symbol_score[domain_idx][symbol]; the
                # `.get()` will default towards the `gap`, which is always necessary
                num_gaps = len([value for value in key if value == self.gap])
                self.scores[key] = np.mean(
                    [
                        symbol_score.get((domain_idx, symbol), num_gaps * gap_score)
                        for domain_idx, symbol in enumerate(key)
                    ]
                )

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
        for sub_key in itertools.product(*sub_alphabets):
            if sub_key not in self.scores:
                if all([value == self.gap for value in sub_key]):
                    continue

                sub_scores = []
                for key, score in self.scores.items():
                    matches = [
                        key[idx] == v for idx, v in enumerate(sub_key) if v is not None
                    ]
                    if all(matches):
                        sub_scores.append(score)

                # Add adjusted value
                # TODO: allow other manipulations? percentile perhaps?
                self.scores[sub_key] = max(sub_scores)

    def _load(self, filename):
        """
        Internal function for loading a serialized matrix from file.

        This function is not supposed to be used directly by the user, but called
        by `__init__()`.
        """

        with open(filename) as json_handler:
            serial_data = json.load(json_handler)

            self.gap = serial_data["gap"]
            self._dr = tuple(serial_data["domain_range"])
            self.alphabets = serial_data["alphabets"].copy()

            # Make sure values are floats
            self.scores = {
                tuple([None if k == "NULL" else k for k in key.split(" / ")]): float(
                    value
                )
                for key, value in serial_data["scores"].items()
            }

    def save(self, filename):
        """
        Serialize a matrix to disk.

        Note that, due limits in serializing with JSON, which do not allow lists
        as keys, the `" / "` symbol (note the spaces) is not allowed.

        Parameters
        ==========
        filename : str
            Path where to write serialized data.
        """

        # Make sure the reserved symbol is not used
        if any([" / " in alphabet for alphabet in self.alphabets]):
            raise ValueError("At least one alphabet uses reserved symbols.")

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

        # Build serialized data
        serial_data = {
            "alphabets": self.alphabets,
            "gap": self.gap,
            "domain_range": list(self._dr),
            "scores": _scores,
        }

        # Open handler and write to disk
        with open(filename, "w") as json_handler:
            json.dump(serial_data, json_handler, indent=4, ensure_ascii=False)

    def compute_submatrices(self, domains):
        """
        Compute submatrices from a collection of domains.

        Sub-matrices are useful, when compared to the possibility of addressing a
        normal matrix with `None` values in keys, in simplifying computation and
        debugging, as well as in some methods for smoothing.
        """

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
        """
        Return a copy of the current matrix object.

        Returns
        =======
        matrix : ScoringMatrix
            A copy of the current ScoringMatrix.
        """
        return copy.deepcopy(self)

    # TODO: currently only working for 2 or 3 domains
    # TODO: use more options from tabulate
    # TODO: check about identity matrix? from alphabets with a method?
    def tabulate(self):
        """
        Build a string with a tabulated representation of the matrix.

        Returns
        =======
        table : str
            A string with the representation of the matrix, intended for human
            consuption.
        """

        rows = []
        if self.domains == 2:

            for symbol_a in self.alphabets[0]:
                row = [symbol_a] + [
                    self.scores[symbol_a, symbol_b] for symbol_b in self.alphabets[1]
                ]
                rows.append(row)

            headers = [""] + list(self.alphabets[1])

        elif self.domains == 3:
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

    def __setitem__(self, key, value):
        # We need to treat `None`, as usual
        matches = [
            k in alphabet for k, alphabet in zip(key, self.alphabets) if k is not None
        ]
        if not all(matches):
            raise ValueError("`key` uses symbol(s) not in alphabet.")

        self.scores[tuple(key)] = value
