# __init__.py

"""
__init__.py file for the `malign` library.
"""

# Package information for the library
__version__ = "0.1.1"
__author__ = "Tiago Tresoldi"
__email__ = "tresoldi@shh.mpg.de"

# Build the namespace
from malign.malign import align
from malign.kbest import fill_scorer, compute_graph, get_aligns

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
