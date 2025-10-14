"""Property-based tests using Hypothesis (Phase 3)."""


# from hypothesis import given, strategies as st
# import malign

# TODO: Phase 3 - Implement property-based tests
# These tests will use hypothesis to generate random inputs and verify properties

# @given(st.lists(st.text(alphabet="ACGT", min_size=1, max_size=20), min_size=2, max_size=5))
# def test_alignment_length_consistency(sequences):
#     """All sequences in alignment should have same length."""
#     alms = malign.align(sequences, k=1)
#     for alm in alms:
#         lengths = [len(seq) for seq in alm.seqs]
#         assert len(set(lengths)) == 1
#
#
# @given(st.lists(st.text(alphabet="ACGT", min_size=1), min_size=2))
# def test_alignment_preserves_symbols(sequences):
#     """Alignment should preserve original symbols (ignoring gaps)."""
#     alms = malign.align(sequences, k=1)
#     for idx, seq in enumerate(sequences):
#         aligned = [s for s in alms[0].seqs[idx] if s != "-"]
#         assert aligned == list(seq)
#
#
# @given(st.lists(st.text(alphabet="ACGT", min_size=1), min_size=2))
# def test_score_ordering(sequences):
#     """Alignments should be ordered by score (best first)."""
#     alms = malign.align(sequences, k=5)
#     scores = [alm.score for alm in alms]
#     assert scores == sorted(scores, reverse=True)
