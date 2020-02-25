#!/usr/bin/env python3

"""
__main__.py

Module for command-line execution of alignment.
"""

# Import Python standard libraries
import argparse

# Import our library
import malign

# TODO: Remove temporary DNA scorer holder in future versions
DNA_SCORER = {
    ("A", "A"): 10,
    ("A", "G"): -1,
    ("A", "C"): -3,
    ("A", "T"): -4,
    ("G", "A"): -1,
    ("G", "G"): 7,
    ("G", "C"): -5,
    ("G", "T"): -3,
    ("C", "A"): -3,
    ("C", "G"): -5,
    ("C", "C"): 9,
    ("C", "T"): 0,
    ("T", "A"): -4,
    ("T", "G"): -3,
    ("T", "C"): 0,
    ("T", "T"): 8,
    ("A", "-"): -5,
    ("G", "-"): -5,
    ("C", "-"): -5,
    ("T", "-"): -5,
    ("-", "A"): -5,
    ("-", "G"): -5,
    ("-", "C"): -5,
    ("-", "T"): -5,
}


def parse_arguments():
    """
    Parses arguments and returns a namespace.
    """

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "seq_a", type=str, help="The first string to be aligned"
    )
    parser.add_argument(
        "seq_b", type=str, help="The second string to be aligned"
    )
    parser.add_argument(
        "--dna",
        action="store_true",
        help="Whether to use the standard DNA scorer (default: uniform scorer)",
    )
    parser.add_argument(
        "-k",
        type=int,
        help="Number of k-best alignments to output (default: 4)",
        default=4,
    )
    args = parser.parse_args()

    return args


def main():
    """
    Entry point for the command-line tool.
    """

    # Parse command-line arguments
    args = parse_arguments()

    # Run alignments
    dna_seq1 = [base for base in args.seq_a]
    dna_seq2 = [base for base in args.seq_b]
    if args.dna:
        scorer = malign.fill_scorer("ACGT", "ACGT", DNA_SCORER)
    else:
        scorer = malign.fill_scorer("ACGT", "ACGT")
    graph = malign.compute_graph(dna_seq1, dna_seq2, scorer)

    dest = "%i:%i" % (len(dna_seq1), len(dna_seq2))
    aligns = malign.get_aligns(graph, ("0:0", dest), dna_seq1, dna_seq2, args.k)

    # Output
    for idx, align in enumerate(aligns):
        print("Alignment #%i (score: %.2f)" % (idx, align[1]))
        print(" ".join(align[0][0]))
        print(" ".join(align[0][1]))
        print()


if __name__ == "__main__":
    main()
