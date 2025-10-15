"""Utilities for loading and working with gold standard test data.

This module provides functions to load cognate sets and alignments from
the YAML files extracted from Arca Verborum gold standard data.
"""

from pathlib import Path
from typing import NamedTuple

import yaml

from malign.alignment import Alignment


class GoldForm(NamedTuple):
    """A single form within a cognate set.

    Attributes:
        language_id: Language identifier.
        glottolog_name: Glottolog name for the language.
        segments: Unaligned segments (original sequence).
        alignment: Aligned segments (with gaps as "-").
    """

    language_id: str
    glottolog_name: str
    segments: list[str]
    alignment: list[str]


class GoldCognateSet(NamedTuple):
    """A cognate set with gold standard alignments.

    Attributes:
        id: Cognate set identifier.
        dataset: Source dataset name.
        gloss: Concepticon gloss (meaning).
        parameter_id: Parameter identifier.
        forms: List of forms in this cognate set.
    """

    id: str
    dataset: str
    gloss: str
    parameter_id: str
    forms: list[GoldForm]


def load_cognate_sets(yaml_path: str | Path) -> list[GoldCognateSet]:
    """Load cognate sets from a YAML file.

    Args:
        yaml_path: Path to YAML file containing cognate sets.

    Returns:
        List of GoldCognateSet objects.

    Example:
        >>> sets = load_cognate_sets("tests/data/cognates/regression_test_set.yml")
        >>> len(sets)
        100
        >>> sets[0].dataset
        'bowernpny'
    """
    with open(yaml_path, encoding="utf-8") as f:
        data = yaml.safe_load(f)

    cognate_sets = []
    for cog_set in data["cognate_sets"]:
        forms = []
        for form in cog_set["forms"]:
            forms.append(
                GoldForm(
                    language_id=form["language_id"],
                    glottolog_name=form["glottolog_name"],
                    segments=form["segments"],
                    alignment=form["alignment"],
                )
            )

        cognate_sets.append(
            GoldCognateSet(
                id=cog_set["id"],
                dataset=cog_set["dataset"],
                gloss=cog_set["gloss"],
                parameter_id=cog_set["parameter_id"],
                forms=forms,
            )
        )

    return cognate_sets


def cognate_set_to_sequences(cog_set: GoldCognateSet) -> list[list[str]]:
    """Convert a cognate set to a list of sequences (segments only).

    Args:
        cog_set: A gold cognate set.

    Returns:
        List of sequences (each sequence is a list of segments).

    Example:
        >>> sets = load_cognate_sets("tests/data/cognates/integration_examples.yml")
        >>> seqs = cognate_set_to_sequences(sets[0])
        >>> len(seqs)
        3
    """
    return [form.segments for form in cog_set.forms]


def cognate_set_to_gold_alignment(cog_set: GoldCognateSet) -> Alignment:
    """Convert a cognate set's gold alignment to an Alignment object.

    Args:
        cog_set: A gold cognate set.

    Returns:
        Alignment object with gold standard aligned sequences.

    Example:
        >>> sets = load_cognate_sets("tests/data/cognates/integration_examples.yml")
        >>> gold_alm = cognate_set_to_gold_alignment(sets[0])
        >>> len(gold_alm.seqs)
        3
    """
    aligned_seqs = [tuple(form.alignment) for form in cog_set.forms]
    return Alignment(aligned_seqs, score=None)  # Gold alignments don't have scores


def filter_cognate_sets(
    cognate_sets: list[GoldCognateSet],
    min_forms: int | None = None,
    max_forms: int | None = None,
    dataset: str | None = None,
    has_gaps: bool | None = None,
) -> list[GoldCognateSet]:
    """Filter cognate sets by various criteria.

    Args:
        cognate_sets: List of cognate sets to filter.
        min_forms: Minimum number of forms (inclusive).
        max_forms: Maximum number of forms (inclusive).
        dataset: Filter by dataset name.
        has_gaps: If True, only sets with gaps; if False, only without gaps.

    Returns:
        Filtered list of cognate sets.

    Example:
        >>> sets = load_cognate_sets("tests/data/cognates/regression_test_set.yml")
        >>> small_sets = filter_cognate_sets(sets, min_forms=2, max_forms=3)
        >>> len(small_sets) < len(sets)
        True
    """
    filtered = cognate_sets

    if min_forms is not None:
        filtered = [cs for cs in filtered if len(cs.forms) >= min_forms]

    if max_forms is not None:
        filtered = [cs for cs in filtered if len(cs.forms) <= max_forms]

    if dataset is not None:
        filtered = [cs for cs in filtered if cs.dataset == dataset]

    if has_gaps is not None:
        if has_gaps:
            # Keep only sets where at least one form has gaps
            filtered = [
                cs for cs in filtered if any("-" in form.alignment for form in cs.forms)
            ]
        else:
            # Keep only sets where no form has gaps
            filtered = [
                cs for cs in filtered if all("-" not in form.alignment for form in cs.forms)
            ]

    return filtered
