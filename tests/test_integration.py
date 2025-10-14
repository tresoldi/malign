"""Integration tests for full alignment pipelines."""

import pytest

import malign


@pytest.mark.integration
def test_full_pipeline_identity_matrix():
    """Test: sequences → identity matrix → align."""
    sequences = ["ACGT", "AGCT", "ATCT"]
    alms = malign.multi_align(sequences, k=2)

    assert len(alms) <= 2
    assert all(len(alm.seqs) == 3 for alm in alms)


@pytest.mark.integration
def test_full_pipeline_custom_matrix():
    """Test: sequences → custom matrix → align."""
    from malign.scoring_matrix import ScoringMatrix

    scores = {
        ("A", "A"): 1.0,
        ("A", "C"): -0.5,
        ("C", "A"): -0.5,
        ("C", "C"): 1.0,
        ("A", "-"): -1.0,
        ("-", "A"): -1.0,
        ("C", "-"): -1.0,
        ("-", "C"): -1.0,
    }
    matrix = ScoringMatrix(scores)

    sequences = ["AC", "AA"]
    alms = malign.multi_align(sequences, k=1, matrix=matrix)
    assert len(alms) == 1


# TODO: Phase 2 - Add more integration tests
# @pytest.mark.integration
# def test_matrix_learning_pipeline():
#     """Test: cognates → learn matrix → align."""
#     pass
#
# @pytest.mark.integration
# def test_yaml_serialization_pipeline():
#     """Test: matrix → save YAML → load → align."""
#     pass
