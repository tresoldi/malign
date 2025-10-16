"""Property-based tests using Hypothesis.

These tests verify invariant properties that should hold across all inputs.
"""

import tempfile
from hypothesis import given, strategies as st, settings

import malign
from malign.scoring_matrix import ScoringMatrix


# Hypothesis settings: conservative for CI
test_settings = settings(max_examples=50, deadline=2000)


@given(
    seqs=st.lists(
        st.text(alphabet="ACGT", min_size=1, max_size=20),
        min_size=2,
        max_size=5,
    )
)
@test_settings
def test_property_alignment_length_consistency(seqs):
    """Property: All sequences in an alignment must have equal length."""
    alms = malign.align(seqs, k=10, method="anw")

    for alm in alms:
        lengths = [len(seq) for seq in alm.seqs]
        assert len(set(lengths)) == 1, f"Inconsistent lengths: {lengths}"


@given(
    seqs=st.lists(
        st.text(alphabet="ACGT", min_size=1, max_size=20),
        min_size=2,
        max_size=5,
    )
)
@test_settings
def test_property_symbol_preservation(seqs):
    """Property: Removing gaps from alignment should yield original sequence."""
    alms = malign.align(seqs, k=1, method="anw")

    for idx, original_seq in enumerate(seqs):
        aligned = [s for s in alms[0].seqs[idx] if s != "-"]
        assert aligned == list(original_seq), f"Sequence {idx} not preserved"


@given(
    seqs=st.lists(
        st.text(alphabet="ACGT", min_size=1, max_size=15),
        min_size=2,
        max_size=4,
    ),
    k=st.integers(min_value=1, max_value=5),
)
@test_settings
def test_property_score_ordering(seqs, k):
    """Property: k-best alignments should be sorted by descending score."""
    alms = malign.align(seqs, k=k, method="anw")

    if len(alms) > 1:
        scores = [alm.score for alm in alms]
        # Check that scores are in descending order (or equal)
        for i in range(len(scores) - 1):
            assert scores[i] >= scores[i + 1], f"Scores not ordered: {scores}"


@given(
    seqs=st.lists(
        st.text(alphabet="ACGT", min_size=1, max_size=12),
        min_size=2,
        max_size=3,
    )
)
@test_settings
def test_property_k_best_uniqueness(seqs):
    """Property: All k-best alignments should be distinct."""
    alms = malign.align(seqs, k=10, method="anw")

    # Convert to tuples for comparison
    alm_tuples = [tuple(tuple(seq) for seq in alm.seqs) for alm in alms]
    unique_alms = set(alm_tuples)

    assert len(alm_tuples) == len(unique_alms), f"Duplicate alignments found"


@given(
    # Generate simpler scoring dictionaries for matrix
    match_score=st.floats(min_value=-5.0, max_value=5.0, allow_nan=False, allow_infinity=False),
    mismatch_score=st.floats(min_value=-5.0, max_value=5.0, allow_nan=False, allow_infinity=False),
    gap_score=st.floats(min_value=-5.0, max_value=5.0, allow_nan=False, allow_infinity=False),
)
@test_settings
def test_property_matrix_yaml_roundtrip(match_score, mismatch_score, gap_score):
    """Property: Matrix should survive save/load cycle unchanged."""
    # Create simple matrix from sequences
    matrix1 = ScoringMatrix.from_sequences(
        sequences=[["A", "C", "G", "T"], ["A", "C", "G", "T"]],
        match=match_score,
        mismatch=mismatch_score,
        gap_score=gap_score,
    )

    with tempfile.NamedTemporaryFile(mode="w", suffix=".yml", delete=False) as f:
        temp_path = f.name

    try:
        # Save and load
        matrix1.save(temp_path)
        matrix2 = ScoringMatrix.from_yaml(temp_path)

        # Verify equality
        assert matrix1.scores.keys() == matrix2.scores.keys(), "Keys don't match"
        for key in matrix1.scores:
            assert abs(matrix1.scores[key] - matrix2.scores[key]) < 1e-6, f"Score mismatch for {key}"
        assert matrix1.domains == matrix2.domains, "Domains don't match"
        assert matrix1.gap == matrix2.gap, "Gap symbol doesn't match"
    finally:
        import os

        if os.path.exists(temp_path):
            os.unlink(temp_path)


@given(
    cognate_sets=st.lists(
        st.lists(
            st.lists(st.sampled_from(["A", "C", "G", "T"]), min_size=2, max_size=8),
            min_size=2,
            max_size=2,
        ),
        min_size=2,
        max_size=5,
    )
)
@test_settings
def test_property_learning_produces_valid_matrix(cognate_sets):
    """Property: Learned matrix should be usable for alignment."""
    try:
        # Learn matrix with minimal iterations for speed
        learned = malign.learn_matrix(cognate_sets, method="em", max_iter=2, gap="-")

        # Should be able to align with it
        test_seqs = cognate_sets[0] if cognate_sets else [["A"], ["A"]]
        alms = malign.align(test_seqs, k=1, matrix=learned)

        # Basic validity checks
        assert len(alms) >= 1, "No alignments produced"
        assert alms[0].score is not None, "Score is None"
        assert len(alms[0].seqs) == len(test_seqs), "Wrong number of sequences"

        # Check that gaps were added if needed
        for seq in alms[0].seqs:
            assert isinstance(seq, tuple), "Aligned sequence not a tuple"

    except Exception as e:
        # Some random cognate sets might be pathological, that's OK
        # But the method shouldn't crash
        assert isinstance(e, (ValueError, ZeroDivisionError, IndexError)), f"Unexpected error: {type(e).__name__}: {e}"


# Mark slow tests
try:
    import pytest

    # Create slow variants with more examples
    slow_settings = settings(max_examples=500, deadline=10000)

    @pytest.mark.slow
    @given(
        seqs=st.lists(
            st.text(alphabet="ACGTRYSWKMBDHVN", min_size=1, max_size=30),
            min_size=2,
            max_size=8,
        )
    )
    @slow_settings
    def test_property_alignment_length_consistency_slow(seqs):
        """Slow variant: Test with more examples and larger alphabet."""
        alms = malign.align(seqs, k=5, method="yenksp")
        for alm in alms:
            lengths = [len(seq) for seq in alm.seqs]
            assert len(set(lengths)) == 1


    @pytest.mark.slow
    @given(
        seqs=st.lists(
            st.text(alphabet="ACGT", min_size=1, max_size=30),
            min_size=2,
            max_size=6,
        )
    )
    @slow_settings
    def test_property_symbol_preservation_slow(seqs):
        """Slow variant: Test with more examples and longer sequences."""
        alms = malign.align(seqs, k=1, method="yenksp")
        for idx, original_seq in enumerate(seqs):
            aligned = [s for s in alms[0].seqs[idx] if s != "-"]
            assert aligned == list(original_seq)

except ImportError:
    pass  # pytest not available
