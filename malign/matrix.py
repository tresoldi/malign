"""
Module for scoring matrices.
"""

# Import Python standard libraries
from collections import defaultdict
import itertools

# Import 3rd-party libraries
import numpy as np

from . import utils

# TODO: add methods for loading and storing the matrices
# TODO: add auxiliary function to build a scoring matrix in a single pass from subs

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
        fill : string or None
            The method used to fill empty elements in the matrix, using
            values computed from the ones provided. If set to `None`, the
            matrix will not be filled. Choices are `"standard"`. Defaults to
            the `"standard"` method.
        """

        # Store values
        self.gap = kwargs.get("gap", "-")

        # Collect submatrices, if they were provided
        sub_matrices = kwargs.get("sub_matrices", {})

        # Make sure all `scores` report the same number of domains, store
        # the number of domains, an set (or override, if necessary) the
        # value of exclusive-gap sites.
        # TODO: check if in line with submatrices
        domains = set([len(key) for key in scores])
        if len(domains) > 1:
            raise ValueError("Different domain-lengths in `scores`.")

        # Copy `domains` and define the global, full `domain range`
        # NOTE: the comma after `self.domains` is used to unpack the
        # `domains` set, which at this point we know contains a single
        # item
        self.domains, = domains
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

        # Extract the alphabets, if they were not provided
        # TODO: check if the alphabets agree with the provided scores
        # TODO: check if the alphabets are list of lists, if provided
        self.alphabets = kwargs.get("alphabets", None)
        if self.alphabets:
            self.alphabets = [sorted(alphabet) for alphabet in self.alphabets]
        else:
            # Collect alphabet from scores
            self.alphabets = [alphabet for alphabet in zip(*scores.keys())]

            # Extend with alphabets from submatrices, if they were provided
            for sub_matrix, obj in sub_matrices.items():
                for dmn, alphabet in zip(sub_matrix, obj.alphabets):
                    self.alphabets[dmn] += alphabet

            # Make sorted lists
            self.alphabets = [
                tuple(sorted(set(alphabet))) for alphabet in self.alphabets
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
        # keys missing indeed; we check for `all` to exclude submatrices
        # that have `None` values
        # TODO: deal with the new Nones (should just filter?), also below
        expected_size = np.prod([len(alphabet) for alphabet in self.alphabets])
        if len([key for key in self.scores if all(key)]) == expected_size:
            return

        # If we have sub-matrices, we will extend them as much as possible
        # to fill the full alignment keys (without overriding any value
        # provided by the user, of course)
        print(self.scores.keys())
        print("REIMPLEMENT WITH NONES")

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
            if key not in self.scores:
                # Collect all sub-scores
                all_sub_scores = []
                for domain_idx in self._dr:
                    sub_key = key[:domain_idx] + key[domain_idx + 1 :]
                    all_sub_scores.append(sub_scores[domain_idx].get(sub_key, None))

                # Cache the new value if possible
                if any(all_sub_scores):
                    score_cache[key] = np.mean([v for v in all_sub_scores if v])

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
                self.scores[key] = np.mean(
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
        sub_alphabets = [
            self.alphabets[d_idx] if d_idx in domain else [None] for d_idx in self._dr
        ]

        # Fill keys
        # TODO: investigate better subsetting, loops to unroll
        # TODO: also what should not be collected like full gaps
        for sub_key in itertools.product(*sub_alphabets):
            sub_scores = []
            for key, score in self.scores.items():
                if all(
                    [key[idx] == v for idx, v in enumerate(sub_key) if v is not None]
                ):
                    sub_scores.append(score)

            # Add adjusted value
            # TODO: allow other manipulations, here just the percentile
            self.scores[sub_key] = np.percentile(sub_scores, 75)

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

        # If there is a domain, compute the new key
        if domain:
            mapper = {d: k for d, k in zip(domain, key)}
            key = tuple([mapper.get(d, None) for d in self._dr])

            # if the key is missing, the submatrix has not been computed yet;
            # we do it now, effectively caching the results
            if key not in self.scores:
                self._fill_domain(domain)

        return self.scores[tuple(key)]
