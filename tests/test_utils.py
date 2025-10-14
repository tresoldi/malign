"""Tests for utility functions."""

import pytest
import malign


def test_tabulate_alms():
    """Test alignment tabulation."""
    alms = malign.multi_align(["ACGT", "AGCT"], k=1)
    table = malign.tabulate_alms(alms)
    assert isinstance(table, str)
    assert len(table) > 0


def test_identity_matrix():
    """Test identity matrix creation."""
    from malign.utils import identity_matrix

    seqs = [["A", "C", "G", "T"], ["A", "C", "G", "T"]]
    matrix = identity_matrix(seqs)
    assert matrix is not None


# TODO: Phase 2/3 - Add tests for new utility functions
# - sum_of_pairs()
# - alignment_entropy()
# - alignment_consensus()
# - alignment_accuracy()
