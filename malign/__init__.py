# __init__.py

"""__init__.py file for the `malign` library."""

# Package information for the library
__version__ = "0.4.0.dev0"
__author__ = "Tiago Tresoldi"
__email__ = "tiago.tresoldi@lingfil.uu.se"

# Build the namespace
from malign.malign import multi_align
from malign.scoring_matrix import ScoringMatrix
from malign.utils import tabulate_alms

# List symbols to export
__all__ = ["ScoringMatrix", "multi_align", "tabulate_alms"]
