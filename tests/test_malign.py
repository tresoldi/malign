#!/usr/bin/env python3

"""
test_malign
===========

Tests for the `malign` package.
"""

# Import Python libraries
import sys
import unittest

# Impor the library itself
import malign


class TestMalign(unittest.TestCase):
    def test_align(self):
        """
        Test pairwise alignment.
        """

        alm_a, alm_b = malign.align("tra", "fata")
        assert tuple(alm_a) == ("f", "a", "t", "a")
        assert tuple(alm_b) == ("t", "r", "a", "-")


if __name__ == "__main__":
    # Explicitly creating and running a test suite allows to profile it
    suite = unittest.TestLoader().loadTestsFromTestCase(TestMalign)
    unittest.TextTestRunner(verbosity=2).run(suite)
