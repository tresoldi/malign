# MAlign API Reference

**Version**: 0.4.0 (in development)

**Note**: This document will be auto-generated from docstrings in Phase 4. Currently serving as a placeholder with structure.

---

## Core Functions

### malign.align()

```python
def align(
    sequences: List[Sequence[Hashable]],
    k: int = 1,
    method: str = "anw",
    matrix: Optional[ScoringMatrix] = None
) -> List[Alignment]
```

Compute multiple alignments for sequences.

**TODO**: Complete documentation in Phase 4 after docstring migration.

---

### malign.align_batch()

```python
def align_batch(
    sequences: List[List[Sequence[Hashable]]],
    k: int = 1,
    method: str = "anw",
    matrix: Optional[ScoringMatrix] = None
) -> List[List[Alignment]]
```

Batch alignment of multiple sequence sets.

**TODO**: Implement in Phase 2, document in Phase 4.

---

## Classes

### Alignment

```python
@dataclass
class Alignment:
    seqs: Sequence[Hashable]
    score: float
```

**TODO**: Add validation methods in Phase 2, document in Phase 4.

---

### ScoringMatrix

```python
class ScoringMatrix:
    def __init__(self, ...):
        ...
```

**TODO**: Refactor with builder pattern in Phase 2, document in Phase 4.

---

## Utility Functions

### malign.tabulate_alms()

```python
def tabulate_alms(alms: List[Alignment]) -> str
```

Return tabulated representation of alignments.

**TODO**: Complete documentation in Phase 4.

---

**Phase 4 Actions**:
1. Auto-generate from docstrings using pydoc-markdown
2. Add examples for each function
3. Add parameter descriptions
4. Add return value descriptions
5. Add cross-references
