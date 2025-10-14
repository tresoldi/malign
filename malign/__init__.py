# __init__.py

"""__init__.py file for the `malign` library."""

# Package information for the library
__version__ = "0.4.0.dev0"
__author__ = "Tiago Tresoldi"
__email__ = "tiago.tresoldi@lingfil.uu.se"

# Build the namespace
from malign.learning import learn_matrix
from malign.malign import multi_align
from malign.metrics import alignment_accuracy, alignment_f1, alignment_precision_recall
from malign.scoring_matrix import ScoringMatrix
from malign.utils import tabulate_alms

# List symbols to export
__all__ = [
    "ScoringMatrix",
    "multi_align",
    "tabulate_alms",
    "alignment_accuracy",
    "alignment_precision_recall",
    "alignment_f1",
    "learn_matrix",
]
