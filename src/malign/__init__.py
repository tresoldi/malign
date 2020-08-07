# __init__.py

"""
__init__.py file for the `malign` library.
"""

# Package information for the library
__version__ = "0.2"
__author__ = "Tiago Tresoldi"
__email__ = "tresoldi@shh.mpg.de"

# Build the namespace
from malign.malign import multi_align

from malign.utils import tabulate_alms
from malign.scoring_matrix import ScoringMatrix
