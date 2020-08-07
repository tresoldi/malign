#!/usr/bin/env python3

"""
__main__.py

Module for command-line execution of alignment.
"""

# Import Python standard libraries
import argparse

# Import our library
import malign


def parse_arguments():
    """
    Parses arguments and returns a namespace.
    """

    parser = argparse.ArgumentParser()
    parser.add_argument("seqs", type=str, help="The sequences to align")
    parser.add_argument(
        "--seq_sep",
        type=str,
        default=",",
        help="The sequence delimiter character (default: `,`)",
    )
    parser.add_argument(
        "--matrix",
        type=str,
        help="Path to the matrix file to use (if not provided, the method default will be used)",
    )
    parser.add_argument(
        "--method",
        type=str,
        default="anw",
        help="The alignment method to use (default: `anw`)",
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

    # Get sequences
    seqs = [[tok for tok in seq] for seq in args.seqs.split(args.seq_sep)]

    # Load matrix, if any
    if args.matrix:
        matrix = malign.ScoringMatrix(filename=args.matrix)
    else:
        matrix = None

    # Run alignment
    alms = malign.multi_align(seqs, args.method, k=args.k, matrix=matrix)

    print(malign.tabulate_alms(alms))


if __name__ == "__main__":
    main()
