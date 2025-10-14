"""Tests for ScoringMatrix builder methods (Phase 2.2)."""

import tempfile
import os

import malign


def test_from_yaml():
    """Test from_yaml() class method."""
    # Create a simple matrix and save it
    scores = {
        ("A", "A"): 1.0,
        ("A", "C"): -0.5,
        ("C", "A"): -0.5,
        ("C", "C"): 1.0,
        ("-", "-"): 0.0,
    }
    matrix1 = malign.ScoringMatrix(scores)

    # Save and load using from_yaml
    with tempfile.NamedTemporaryFile(mode="w", suffix=".yml", delete=False) as f:
        temp_path = f.name

    try:
        matrix1.save(temp_path)
        matrix2 = malign.ScoringMatrix.from_yaml(temp_path)

        # Verify they match
        assert matrix1.scores == matrix2.scores
        assert matrix1.domains == matrix2.domains
        assert matrix1.gap == matrix2.gap
    finally:
        if os.path.exists(temp_path):
            os.unlink(temp_path)


def test_from_sequences_basic():
    """Test from_sequences() with basic parameters."""
    matrix = malign.ScoringMatrix.from_sequences(
        sequences=[["A", "C"], ["X", "Y"]],
        match=1.0,
        mismatch=-0.5,
        gap_score=-1.0,
    )

    # Check structure
    assert matrix.num_domains == 2
    assert matrix.gap == "-"
    assert len(matrix.domains) == 2

    # Check scores
    assert matrix["-", "-"] == 0.0  # All-gap vector
    assert matrix["A", "-"] == -1.0  # Gap score
    assert matrix["-", "X"] == -1.0  # Gap score
    assert matrix["A", "Y"] == -0.5  # Mismatch (different symbols)
    # No natural matches since alphabets are disjoint


def test_from_sequences_custom_gap():
    """Test from_sequences() with custom gap symbol."""
    matrix = malign.ScoringMatrix.from_sequences(
        sequences=[["A", "T"], ["G", "C"]],
        match=2.0,
        mismatch=-1.0,
        gap=".",
        gap_score=-2.0,
    )

    assert matrix.gap == "."
    assert matrix[".", "."] == 0.0
    assert matrix["A", "."] == -2.0


def test_from_sequences_identity():
    """Test from_sequences() for identity-style matrix (same alphabet)."""
    # When using same symbols in both domains, we can get match scores
    matrix = malign.ScoringMatrix.from_sequences(
        sequences=[["A", "C", "G", "T"], ["A", "C", "G", "T"]],
        match=1.0,
        mismatch=-0.5,
        gap_score=-1.0,
    )

    # Same symbol in both domains = match
    assert matrix["A", "A"] == 1.0
    assert matrix["C", "C"] == 1.0

    # Different symbols = mismatch
    assert matrix["A", "C"] == -0.5
    assert matrix["C", "A"] == -0.5
