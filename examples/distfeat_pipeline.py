"""Full pipeline example: distfeat features -> scoring matrix -> alignment -> evaluation.

This example demonstrates the JOSS paper story: using phonological feature
distances from distfeat to build a linguistically-informed scoring matrix,
then aligning cognate sequences and evaluating the results.

Requires: pip install malign[features]
"""

import malign


def main():
    # Step 1: Define cognate sequences (Italian vs Spanish cognates)
    # These are IPA-like phonological segments
    cognate_pairs = [
        (["n", "o", "t", "e"], ["n", "o", "tʃ", "e"]),  # notte / noche
        (["f", "a", "t", "o"], ["h", "a", "d", "o"]),  # fatto / hado
        (["p", "j", "e", "n", "o"], ["ʎ", "e", "n", "o"]),  # pieno / lleno
    ]

    # Step 2: Build scoring matrix from phonological feature distances
    # Collect all unique segments per position
    seq_a_segments = sorted({s for pair in cognate_pairs for s in pair[0]})
    seq_b_segments = sorted({s for pair in cognate_pairs for s in pair[1]})

    matrix = malign.ScoringMatrix.from_distfeat(
        sequences=[seq_a_segments, seq_b_segments],
        gap="-",
        gap_score=-1.0,
    )

    print("Scoring matrix (feature-based):")
    print(matrix.tabulate())
    print()

    # Step 3: Align each cognate pair
    for seq_a, seq_b in cognate_pairs:
        alms = malign.align([seq_a, seq_b], k=3, matrix=matrix, method="anw")

        print(f"Aligning: {' '.join(seq_a)} ~ {' '.join(seq_b)}")
        print(malign.tabulate_alms(alms[:2]))
        print()

    # Step 4: Evaluate alignment quality against gold standard
    gold = malign.alignment.Alignment(
        [("n", "o", "t", "e"), ("n", "o", "tʃ", "e")],
        score=0.0,
    )
    predicted = malign.align(
        [["n", "o", "t", "e"], ["n", "o", "tʃ", "e"]],
        k=1,
        matrix=matrix,
    )[0]

    acc = malign.alignment_accuracy(predicted, gold)
    prec, rec = malign.alignment_precision_recall(predicted, gold)
    f1 = malign.alignment_f1(predicted, gold)

    print(f"Accuracy:  {acc:.2%}")
    print(f"Precision: {prec:.2%}")
    print(f"Recall:    {rec:.2%}")
    print(f"F1 Score:  {f1:.2%}")


if __name__ == "__main__":
    main()
