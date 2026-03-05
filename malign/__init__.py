"""malign: Multiple asymmetric sequence alignment for computational linguistics."""

__version__ = "0.5.0"
__author__ = "Tiago Tresoldi"
__email__ = "tiago.tresoldi@lingfil.uu.se"

from malign.alignment import Alignment
from malign.learning import learn_matrix
from malign.malign import align
from malign.metrics import alignment_accuracy, alignment_f1, alignment_precision_recall
from malign.scoring_matrix import ScoringMatrix
from malign.utils import tabulate_alms

__all__ = [
    "Alignment",
    "ScoringMatrix",
    "align",
    "alignment_accuracy",
    "alignment_f1",
    "alignment_precision_recall",
    "learn_matrix",
    "tabulate_alms",
]
