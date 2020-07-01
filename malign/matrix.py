"""
Module for scoring matrices.
"""

# Import Python standard libraries
from collections import defaultdict
import itertools

# Import 3rd-party libraries
import numpy as np

import utils


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
        self.scores = scores
        self.gap = kwargs.get("gap", "-")

        # Make sure all `scores` report the same number of domains, store
        # the number of domains, an set (or override, if necessary) the
        # value of exclusive-gap sites.
        domains = set([len(key) for key in self.scores])
        if len(domains) > 1:
            raise ValueError("Different domain-lengths in `scores`.")
        else:
            # NOTE: the comma after `self.domains` is used to unpack the
            # `domains` set, which at this point we know contains a single item
            self.domains, = domains
            self.scores[tuple([self.gap] * self.domains)] = 0.0

        # Extract the alphabets, if they were not provided
        # TODO: check if the alphabets agree with the provided scores
        # TODO: check if the alphabets are list of lists, if provided
        self.alphabets = kwargs.get("alphabets", None)
        if self.alphabets:
            self.alphabets = [sorted(alphabet) for alphabet in self.alphabets]
        else:
            self.alphabets = [sorted(set(alphabet)) for alphabet in zip(*scores.keys())]

        # Fill the matrix with the appropriate method if requested
        fill = kwargs.get("fill", "standard")
        if fill:
            if fill == "standard":
                self._fill_matrix(fill)
            else:
                raise ValueError("Unknown filling method.")

    # TODO: currently disregarding the `method`, as there is a single one
    # TODO: replace the `enumerate(self.alphabets)` with a cached range
    # TODO: describe how the standard method can be considered a kind of MLE
    def _fill_matrix(self, method):
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
        if len(self.scores) == expected_size:
            return

        # For each domain, we first collect all keys without considering
        # the symbol in the domain itself, and fill any missing spot with
        # the mean value. No difference is made in terms of gaps.
        # Example: if we have ('a', 'A', '1')=10 and ('b', 'A', '1')=0, but
        # no ('c', 'A', '1'), this will set the latter to the mean value
        # of 5, as coming from (?, 'A', '1')
        sub_scores = {}
        for domain_idx, alphabet in enumerate(self.alphabets):
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
                for domain_idx, alphabet in enumerate(self.alphabets):
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
        if len(self.scores) == expected_size:
            return

        symbol_score = defaultdict(lambda: defaultdict(list))
        for key, score in self.scores.items():
            for domain_idx, alphabet in enumerate(self.alphabets):
                symbol_score[domain_idx][key[domain_idx]].append(score)

        for domain_idx, alphabet in enumerate(self.alphabets):
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

    def __getitem__(self, key):
        """
        Return the score associated with a tuple of alignments.
        
        Parameters
        ==========

        key : tuple of strings
            The element in the multidimensional matrix to be queried.

        Returns
        =======

        value : float
            The value associated with the element.
        """

        return self.scores[key]


if __name__ == "__main__":
    m = ScoringMatrix(utils.DNA_MATRIX)

    s = {
        ("a", "A", "1"): -1,
        ("a", "A", "2"): 4,
        ("a", "A", "3"): 3,
        ("a", "B", "1"): 1,
        ("b", "A", "1"): -10,
        ("b", "A", "2"): 10,
        ("b", "A", "3"): 10,
        ("b", "A", "4"): 10,
        ("c", "B", "1"): 2,
        ("c", "B", "2"): 2,
        ("a", "-", "-"): -2,
        ("b", "-", "-"): -2,
        ("-", "A", "-"): -20,
        ("-", "B", "-"): -20,
        ("-", "-", "1"): -3,
        ("a", "B", "-"): -10,
        ("-", "A", "1"): -100,
        ("-", "A", "2"): -10,
        ("-", "B", "3"): -5,
    }
    m = ScoringMatrix(s)
