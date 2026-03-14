"""Shared pytest fixtures for the malign test suite."""

from pathlib import Path

import pytest

import malign


# --- Cognate data fixtures ---


@pytest.fixture
def dna_cognates():
    """Synthetic DNA cognate sets for testing."""
    return [
        [["A", "C", "G", "T"], ["A", "C", "G", "T"]],
        [["A", "C", "C", "T"], ["A", "G", "C", "T"]],
        [["T", "G", "A", "C"], ["T", "G", "A", "C"]],
        [["G", "A", "T", "T"], ["G", "A", "A", "T"]],
    ]


@pytest.fixture
def identity_matrix_fixture():
    """Factory fixture for creating identity-style scoring matrices."""

    def _make(seqs, **kwargs):
        kwargs.setdefault("match", 1.0)
        kwargs.setdefault("gap_score", -1.0)
        return malign.ScoringMatrix.from_sequences(seqs, **kwargs)

    return _make


# --- Path fixtures ---


@pytest.fixture
def data_dir():
    """Path to the test data/cognates directory."""
    return Path(__file__).parent / "data" / "cognates"


@pytest.fixture
def regression_set(data_dir):
    """Path to the regression test set YAML."""
    return data_dir / "regression_test_set.yml"


@pytest.fixture
def training_set(data_dir):
    """Path to the learning training set YAML."""
    return data_dir / "learning_training_set.yml"


@pytest.fixture
def eval_set(data_dir):
    """Path to the learning evaluation set YAML."""
    return data_dir / "learning_eval_set.yml"
