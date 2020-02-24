#!/usr/bin/env python3

"""
test_malign
===========

Tests for the `malign` package.
"""

# Import Python libraries
import math
import sys
import unittest

# Impor the library itself
import malign


class TestMalign(unittest.TestCase):
    def test_dumb_align(self):
        """
        Test `dumb` pairwise alignment.
        """

        alm_a, alm_b, score = malign.align("tra", "fata", method="dumb")
        assert tuple(alm_a) == ("t", "r", "a", "-")
        assert tuple(alm_b) == ("f", "a", "t", "a")
        assert math.isclose(score, 0.75)

    def test_nw_align(self):
        """
        Test `nw` pairwise alignment.
        """

        alm_a, alm_b, score = malign.align("tra", "fata", method="nw")
        assert tuple(alm_a) == ("-", "-", "t", "r", "a")
        assert tuple(alm_b) == ("f", "a", "t", "-", "a")
        assert math.isclose(score, 0.0)


if __name__ == "__main__":
    # Explicitly creating and running a test suite allows to profile it
    suite = unittest.TestLoader().loadTestsFromTestCase(TestMalign)
    unittest.TextTestRunner(verbosity=2).run(suite)
