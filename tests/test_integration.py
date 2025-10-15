"""Integration tests for full alignment pipelines."""

import tempfile
from pathlib import Path

import pytest

import malign
from malign.alignment import Alignment
from malign.metrics import alignment_accuracy, alignment_f1, alignment_precision_recall
from malign.scoring_matrix import ScoringMatrix


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


@pytest.mark.integration
def test_learning_pipeline_end_to_end():
    """Integration 1: cognate sets → learn matrix → align → validate.

    Tests complete learning workflow from cognate sets to final alignment,
    verifying that learned matrices produce valid alignments.
    """
    # Prepare cognate sets for learning
    cognate_sets = [
        [["A", "C", "G"], ["A", "C", "G"]],  # Perfect match
        [["A", "C", "G", "T"], ["A", "G", "G", "T"]],  # One mismatch
        [["T", "G", "C"], ["T", "G", "C"]],  # Perfect match
        [["C", "C", "A"], ["C", "G", "A"]],  # One mismatch
    ]

    # Learn matrix from cognate sets
    learned_matrix = malign.learn_matrix(cognate_sets, method="em", max_iter=5, gap="-")

    # Verify learned matrix is valid
    assert isinstance(learned_matrix, ScoringMatrix)
    assert len(learned_matrix.scores) > 0
    assert learned_matrix.gap == "-"

    # Use learned matrix to align new sequences
    test_sequences = [["A", "C", "G", "T"], ["A", "C", "G", "T"]]
    alms = malign.multi_align(test_sequences, k=3, matrix=learned_matrix, method="anw")

    # Verify alignments are valid
    assert len(alms) >= 1
    assert len(alms[0].seqs) == 2
    assert alms[0].score is not None

    # Verify all sequences have same length
    for alm in alms:
        lengths = [len(seq) for seq in alm.seqs]
        assert len(set(lengths)) == 1

    # Verify original sequences are preserved (gaps removed)
    for alm in alms:
        for i, original in enumerate(test_sequences):
            aligned_no_gaps = [s for s in alm.seqs[i] if s != "-"]
            assert aligned_no_gaps == original


@pytest.mark.integration
def test_matrix_persistence_pipeline():
    """Integration 2: create → save YAML → load → align → verify identical results.

    Tests matrix serialization round-trip ensuring saved/loaded matrices
    produce identical alignments.
    """
    # Create initial matrix
    scores = {
        ("A", "A"): 2.0,
        ("A", "C"): -1.0,
        ("C", "A"): -1.0,
        ("C", "C"): 2.0,
        ("A", "G"): -0.5,
        ("G", "A"): -0.5,
        ("G", "G"): 2.0,
        ("C", "G"): -0.5,
        ("G", "C"): -0.5,
        ("A", "-"): -2.0,
        ("-", "A"): -2.0,
        ("C", "-"): -2.0,
        ("-", "C"): -2.0,
        ("G", "-"): -2.0,
        ("-", "G"): -2.0,
    }
    original_matrix = ScoringMatrix(scores, gap="-")

    # Align with original matrix
    test_sequences = [["A", "C", "G"], ["A", "G", "G"]]
    alms_original = malign.multi_align(test_sequences, k=5, matrix=original_matrix, method="anw")

    # Save matrix to YAML
    with tempfile.NamedTemporaryFile(mode="w", suffix=".yml", delete=False) as f:
        temp_path = f.name

    try:
        original_matrix.save(temp_path)

        # Load matrix from YAML
        loaded_matrix = ScoringMatrix.from_yaml(temp_path)

        # Verify matrix integrity
        assert loaded_matrix.gap == original_matrix.gap
        assert loaded_matrix.domains == original_matrix.domains
        assert loaded_matrix.scores.keys() == original_matrix.scores.keys()
        for key in original_matrix.scores:
            assert abs(loaded_matrix.scores[key] - original_matrix.scores[key]) < 1e-10

        # Align with loaded matrix
        alms_loaded = malign.multi_align(test_sequences, k=5, matrix=loaded_matrix, method="anw")

        # Verify identical results
        assert len(alms_original) == len(alms_loaded)
        for orig_alm, loaded_alm in zip(alms_original, alms_loaded):
            assert orig_alm.seqs == loaded_alm.seqs
            assert abs(orig_alm.score - loaded_alm.score) < 1e-10

    finally:
        Path(temp_path).unlink(missing_ok=True)


@pytest.mark.integration
def test_anw_vs_yenksp_consistency():
    """Integration 3: ANW vs YenKSP consistency.

    Tests that both alignment methods produce consistent top-k results,
    verifying algorithm correctness across implementations.
    """
    # Test with sequences where both methods should agree on best alignment
    sequences = [["A", "C", "G", "T"], ["A", "G", "G", "T"], ["A", "C", "C", "T"]]

    # Get top-k alignments from both methods
    k = 3
    alms_anw = malign.multi_align(sequences, k=k, method="anw")
    alms_yenksp = malign.multi_align(sequences, k=k, method="yenksp")

    # Both should produce results
    assert len(alms_anw) >= 1
    assert len(alms_yenksp) >= 1

    # Both should have same number of sequences
    assert len(alms_anw[0].seqs) == len(alms_yenksp[0].seqs) == 3

    # Top alignment should be very similar (within tolerance for floating point)
    # Best alignment from both should have similar scores
    score_diff = abs(alms_anw[0].score - alms_yenksp[0].score)
    assert score_diff < 0.5, f"Top scores differ too much: ANW={alms_anw[0].score}, YenKSP={alms_yenksp[0].score}"

    # All alignments should preserve original sequences
    for alms, method_name in [(alms_anw, "ANW"), (alms_yenksp, "YenKSP")]:
        for alm in alms:
            for i, original in enumerate(sequences):
                aligned_no_gaps = [s for s in alm.seqs[i] if s != "-"]
                assert aligned_no_gaps == original, f"{method_name} failed to preserve sequence {i}"

            # All sequences should have same length
            lengths = [len(seq) for seq in alm.seqs]
            assert len(set(lengths)) == 1, f"{method_name} produced inconsistent lengths"


@pytest.mark.integration
def test_metrics_validation_chain():
    """Integration 4: align → compute all metrics → verify mathematical relationships.

    Tests alignment quality metrics, verifying mathematical invariants:
    - Accuracy ∈ [0, 1]
    - Precision, Recall ∈ [0, 1]
    - F1 = 2·P·R/(P+R)
    - Perfect match: all metrics = 1.0
    """
    # Create predicted and gold alignments
    # Case 1: Perfect match
    perfect_pred = Alignment(
        [("A", "C", "G", "T"), ("A", "C", "G", "T"), ("A", "C", "G", "T")], score=10.0
    )
    perfect_gold = Alignment(
        [("A", "C", "G", "T"), ("A", "C", "G", "T"), ("A", "C", "G", "T")], score=10.0
    )

    acc_perfect = alignment_accuracy(perfect_pred, perfect_gold)
    prec_perfect, rec_perfect = alignment_precision_recall(perfect_pred, perfect_gold)
    f1_perfect = alignment_f1(perfect_pred, perfect_gold)

    # Perfect match: all should be 1.0
    assert acc_perfect == 1.0
    assert prec_perfect == 1.0
    assert rec_perfect == 1.0
    assert f1_perfect == 1.0

    # Case 2: Partial match
    partial_pred = Alignment(
        [("A", "C", "-", "G", "T"), ("A", "C", "G", "-", "T"), ("A", "C", "G", "G", "T")],
        score=5.0,
    )
    partial_gold = Alignment(
        [("A", "C", "G", "-", "T"), ("A", "C", "-", "G", "T"), ("A", "C", "G", "G", "T")],
        score=5.0,
    )

    acc_partial = alignment_accuracy(partial_pred, partial_gold)
    prec_partial, rec_partial = alignment_precision_recall(partial_pred, partial_gold)
    f1_partial = alignment_f1(partial_pred, partial_gold)

    # All metrics should be in [0, 1]
    assert 0.0 <= acc_partial <= 1.0
    assert 0.0 <= prec_partial <= 1.0
    assert 0.0 <= rec_partial <= 1.0
    assert 0.0 <= f1_partial <= 1.0

    # Verify F1 = 2·P·R/(P+R)
    if prec_partial + rec_partial > 0:
        expected_f1 = 2 * prec_partial * rec_partial / (prec_partial + rec_partial)
        assert abs(f1_partial - expected_f1) < 1e-10

    # Partial match should have lower metrics than perfect
    assert acc_partial < acc_perfect
    assert f1_partial < f1_perfect

    # Case 3: Partial mismatch (columns 0 and 3 match, 1 and 2 don't)
    mismatch_pred = Alignment([("A", "C", "G", "T"), ("A", "G", "C", "T")], score=1.0)
    mismatch_gold = Alignment([("A", "T", "T", "T"), ("A", "T", "T", "T")], score=1.0)

    acc_mismatch = alignment_accuracy(mismatch_pred, mismatch_gold)

    # Column 0: (A,A) vs (A,A) = match
    # Column 1: (C,G) vs (T,T) = no match
    # Column 2: (G,C) vs (T,T) = no match
    # Column 3: (T,T) vs (T,T) = match
    # Result: 2/4 = 0.5
    assert acc_mismatch == 0.5  # 2 out of 4 columns match
    assert acc_mismatch < acc_perfect


@pytest.mark.slow
@pytest.mark.integration
def test_large_scale_stress_test():
    """Integration 5: Large-scale stress test.

    Tests system behavior with many sequences and multiple cognate sets,
    verifying robustness under realistic workload.

    Marked @pytest.mark.slow as this can take several seconds.
    """
    # Test 1: Many sequences (stress alignment algorithm)
    large_sequence_set = [
        ["A", "C", "G", "T"],
        ["A", "C", "G", "T"],
        ["A", "G", "G", "T"],
        ["A", "C", "C", "T"],
        ["T", "C", "G", "T"],
        ["A", "C", "G", "G"],
    ]

    # Should handle 6 sequences (reduced for reasonable runtime)
    alms_large = malign.multi_align(large_sequence_set, k=2, method="anw")
    assert len(alms_large) >= 1
    assert len(alms_large[0].seqs) == 6

    # Verify consistency
    for alm in alms_large:
        lengths = [len(seq) for seq in alm.seqs]
        assert len(set(lengths)) == 1

    # Test 2: Many cognate sets (stress learning algorithm)
    many_cognate_sets = []
    for i in range(20):  # 20 cognate sets (reduced for reasonable runtime)
        # Generate simple cognate sets with variation
        base = ["A", "C", "G", "T"]
        variant = base.copy()
        variant[i % 4] = "G"  # Introduce variation
        many_cognate_sets.append([base, variant])

    # Should handle learning from many sets
    learned = malign.learn_matrix(many_cognate_sets, method="em", max_iter=2, gap="-")
    assert isinstance(learned, ScoringMatrix)
    assert len(learned.scores) > 0

    # Test 3: High k value (stress k-best algorithm)
    test_seqs = [["A", "C", "G"], ["A", "G", "G"], ["T", "C", "G"]]
    alms_high_k = malign.multi_align(test_seqs, k=20, method="anw")

    # Should return up to k alignments (may be fewer if not enough exist)
    assert 1 <= len(alms_high_k) <= 20

    # All alignments should be distinct
    alm_tuples = [tuple(tuple(seq) for seq in alm.seqs) for alm in alms_high_k]
    assert len(alm_tuples) == len(set(alm_tuples)), "Duplicate alignments found"

    # Scores should be in descending order
    scores = [alm.score for alm in alms_high_k]
    for i in range(len(scores) - 1):
        assert scores[i] >= scores[i + 1], f"Scores not ordered: {scores}"

    # Test 4: Long sequences (stress memory and computation)
    long_seq1 = ["A", "C", "G", "T"] * 6  # 24 symbols (reduced for reasonable runtime)
    long_seq2 = ["A", "C", "G", "T"] * 6
    long_seq2[5] = "G"  # Introduce variation

    alms_long = malign.multi_align([long_seq1, long_seq2], k=3, method="anw")
    assert len(alms_long) >= 1

    # Verify sequences preserved
    for i, original in enumerate([long_seq1, long_seq2]):
        aligned_no_gaps = [s for s in alms_long[0].seqs[i] if s != "-"]
        assert aligned_no_gaps == original
