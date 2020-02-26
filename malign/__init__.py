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
from malign.nw import *
from malign.kbest import fill_scorer, compute_graph, get_aligns

from malign.utils import DNA_SCORER, print_alms
