# __init__.py

"""__init__.py file for the `malign` library."""

# Package information for the library
__version__ = "0.4.0b1"
__author__ = "Tiago Tresoldi"
__email__ = "tiago.tresoldi@lingfil.uu.se"

# Build the namespace
from malign.learning import learn_matrix
from malign.malign import align
from malign.metrics import alignment_accuracy, alignment_f1, alignment_precision_recall
from malign.scoring_matrix import ScoringMatrix
from malign.utils import tabulate_alms

# List symbols to export
__all__ = [
    "ScoringMatrix",
    "align",
    "alignment_accuracy",
    "alignment_f1",
    "alignment_precision_recall",
    "learn_matrix",
    "tabulate_alms",
]
