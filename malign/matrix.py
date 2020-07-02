"""
Module for scoring matrices.
"""

# Import Python standard libraries
from collections import defaultdict
import itertools

# Import 3rd-party libraries
import numpy as np

from . import utils

# TODO: add functions for loading and storing the matrices


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

    def __init__(self, scores, **kwargs):
        """
        Initialize a scoring matrix.

        Parameters
        ==========

        scores : dict
            A scoring dictionary, with tuples of strings as keys (respecting
            the order that will be used when using the matrix) and floating
            points as values. Exclusive-gap alignments (that is, keys
            composed only of gaps) will be overriden to a value of 0.0 if
            provided.
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
        fill : string or None
            The method used to fill empty elements in the matrix, using
            values computed from the ones provided. If set to `None`, the
            matrix will not be filled. Choices are `"standard"`. Defaults to
            the `"standard"` method.
        """

        # Store values
        self.gap = kwargs.get("gap", "-")

        # If `scores` is a dictionary, it is dictionary (potentially
        # incomplete and which needs to be filled) of full alignment
        # sites
        if isinstance(scores, dict):
            # Make sure all `scores` report the same number of domains, store
            # the number of domains, an set (or override, if necessary) the
            # value of exclusive-gap sites.
            domains = set([len(key) for key in scores])
            if len(domains) > 1:
                raise ValueError("Different domain-lengths in `scores`.")

            # Copy `domains` and define the global, full `domain range`
            # NOTE: the comma after `self.domains` is used to unpack the
            # `domains` set, which at this point we know contains a single
            # item
            self.domains, = domains
            self._dr = tuple(range(self.domains))
            self.scores = {self._dr: scores.copy()}
            self.scores[self._dr][tuple([self.gap] * self.domains)] = 0.0

            # Extract the alphabets, if they were not provided
            # TODO: check if the alphabets agree with the provided scores
            # TODO: check if the alphabets are list of lists, if provided
            self.alphabets = kwargs.get("alphabets", None)
            if self.alphabets:
                self.alphabets = [sorted(alphabet) for alphabet in self.alphabets]
            else:
                self.alphabets = [
                    sorted(set(alphabet)) for alphabet in zip(*scores.keys())
                ]

            # Fill the matrix with the appropriate method if requested
            fill = kwargs.get("fill", "standard")
            if fill:
                if fill == "standard":
                    self._fill_full_matrix(fill)
                else:
                    raise ValueError("Unknown filling method.")

    # TODO: currently disregarding the `method`, as there is a single one
    # TODO: describe how the standard method can be considered a kind of MLE
    # TODO: use `self.domain_range`
    def _fill_full_matrix(self, method):
        """
        Internal function for filling a matrix if there are missing values.

        Parameters
        ==========

        method : string
            The method to be used for filling the matrix. Choices are
            `"standard"`.
        """

        # Compute the expected size of the scorer, so we know if there are
        # keys missing indeed
        expected_size = np.prod([len(alphabet) for alphabet in self.alphabets])
        if len(self.scores[self._dr]) == expected_size:
            return

        # Cache the range of domains
        domains = list(range(len(self.alphabets)))

        # For each domain, we first collect all keys without considering
        # the symbol in the domain itself, and fill any missing spot with
        # the mean value. No difference is made in terms of gaps.
        # Example: if we have ('a', 'A', '1')=10 and ('b', 'A', '1')=0, but
        # no ('c', 'A', '1'), this will set the latter to the mean value
        # of 5, as coming from (?, 'A', '1')
        sub_scores = {}
        for domain_idx in domains:
            # Collect all sub-keys and values
            sub_scores[domain_idx] = defaultdict(list)
            for key, score in self.scores[self._dr].items():
                sub_key = tuple(
                    [value for idx, value in enumerate(key) if idx != domain_idx]
                )
                sub_scores[domain_idx][sub_key].append(score)

        # Take the mean value of all sub_scores
        for domain_idx, scores in sub_scores.items():
            sub_scores[domain_idx] = {
                sub_key: np.mean(values) for sub_key, values in scores.items()
            }

        # Set the new score when possible, from the mean of the
        # available `sub_scores`. We cache the new values in order to only
        # apply them after the loop, not changing `self.scores` in-place.
        score_cache = {}
        for key in itertools.product(*self.alphabets):
            if key not in self.scores[self._dr]:
                # Collect all sub-scores
                all_sub_scores = []
                for domain_idx in domains:
                    sub_key = key[:domain_idx] + key[domain_idx + 1 :]
                    all_sub_scores.append(sub_scores[domain_idx].get(sub_key, None))

                # Cache the new value if possible
                if any(all_sub_scores):
                    score_cache[key] = np.mean([v for v in all_sub_scores if v])

        # Update with the new values
        self.scores[self._dr].update(score_cache)

        # If there are still keys missing, collect the mean sub_score
        # by symbol in each domain and set to this value -- this is similar
        # to what we do above, but much less precise as we take all symbols
        # independently
        if len(self.scores[self._dr]) == expected_size:
            return

        symbol_score = defaultdict(lambda: defaultdict(list))
        for key, score in self.scores[self._dr].items():
            for domain_idx in domains:
                symbol_score[domain_idx][key[domain_idx]].append(score)

        for domain_idx in domains:
            symbol_score[domain_idx] = {
                symbol: np.mean(scores)
                for symbol, scores in symbol_score[domain_idx].items()
            }

        for key in itertools.product(*self.alphabets):
            if key not in self.scores[self._dr]:
                print("!!!!", key)
                self.scores[self._dr][key] = np.mean(
                    [
                        symbol_score[domain_idx][symbol]
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
        sub_alphabets = [self.alphabets[d_idx] for d_idx in domain]

        # Fill keys
        # TODO: investigate better subsetting, loops to unroll
        # TODO: also what should not be collected like full gaps
        self.scores[domain] = {}
        for sub_key in itertools.product(*sub_alphabets):
            # Collect all values in the main matrix that match the sub_key
            idx_values = list(zip(domain, sub_key))
            sub_scores = []
            for key, score in self.scores[self._dr].items():
                if all([key[idx] == v for idx, v in idx_values]):
                    sub_scores.append(score)

            # Add adjusted value
            # TODO: allow other manipulations, here just the percentile
            self.scores[domain][sub_key] = np.percentile(sub_scores, 75)

        # Fill full gap
        self.scores[domain][tuple(self.gap * len(domain))] = 0.0

    def __call__(self, key, domain=None):
        """
        Return the score associated with a tuple of alignments per domain.

        If a sub-domain is requested and it has not been computed so far,
        it will be inferred and cached for reuse.

        Parameters
        ==========

        key : tuple of strings
            The element in the multidimensional matrix to be queried.
        domain : tuple or list of integers
            The domain where the key should be queried. Defaults to the
            complete alignment site.

        Returns
        =======

        value : float
            The value associated with the element.
        """

        if not domain:
            return self.scores[self._dr][tuple(key)]

        # Infer and cache the domain if not available
        domain = tuple(domain)
        if domain not in self.scores:
            self._fill_domain(domain)

        return self.scores[domain][tuple(key)]
