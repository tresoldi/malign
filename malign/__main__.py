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
    parser.add_argument(
        "seq_a", type=str, help="The first string to be aligned"
    )
    parser.add_argument(
        "seq_b", type=str, help="The second string to be aligned"
    )
    args = parser.parse_args()

    return args


def main():
    """
    Entry point for the command-line tool.
    """

    # Parse command-line arguments
    args = parse_arguments()

    # Get alignments
    alm_a, alm_b = malign.align(args.seq_a, args.seq_b)

    # Output alignments
    print("A: [%s]" % " ".join(alm_a))
    print("B: [%s]" % " ".join(alm_b))


if __name__ == "__main__":
    main()
