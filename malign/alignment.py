"""Module for the Alignment class.

The `Alignment` class is a simple data class that holds aligned sequences and
their score. It was originally a dictionary passed back and forth among
functions, for which a data class is a good replacement.
"""

from collections.abc import Hashable, Sequence
from dataclasses import dataclass

# TODO: write methods for comparison, based on score
# TODO: add various checks post-initialization


@dataclass
class Alignment:
    """Data class for holding aligned sequences and their score.

    @param seqs: The sequences in the alignment.
    @param score: The alignment score.
    """

    seqs: Sequence[Sequence[Hashable]]
    score: float

    def __len__(self) -> int:
        """Return the number of sequences in the alignment.

        @return: The number of sequences in the alignment.
        """
        return len(self.seqs)

    def __getitem__(self, idx: int) -> Sequence[Hashable]:
        """Return a sequence by its index.

        @param idx: The index of the sequence in the alignment.
        @return: The sequence at the requested index.
        """
        return self.seqs[idx]
