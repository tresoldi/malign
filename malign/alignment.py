"""Module for the Alignment dataclass."""

from collections.abc import Hashable, Sequence
from dataclasses import dataclass


@dataclass
class Alignment:
    """Data class for holding aligned sequences and their score.

    Attributes:
        seqs: The sequences in the alignment.
        score: The alignment score.
    """

    seqs: Sequence[Sequence[Hashable]]
    score: float

    def __len__(self) -> int:
        """Return the number of sequences in the alignment."""
        return len(self.seqs)

    def __getitem__(self, idx: int) -> Sequence[Hashable]:
        """Return a sequence by its index."""
        return self.seqs[idx]
