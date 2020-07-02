# __init__.py

"""
__init__.py file for the `malign` library.
"""

# Package information for the library
__version__ = "0.1.1"
__author__ = "Tiago Tresoldi"
__email__ = "tresoldi@shh.mpg.de"

# Build the namespace
# TODO: be consistent in importing
from malign.malign import pw_align, multi_align

# from malign.nw import *
# from malign.graph import *

from malign.utils import print_alms, print_malms
from malign.matrix import ScoringMatrix
