# MAlign

[![PyPI](https://img.shields.io/pypi/v/malign.svg)](https://pypi.org/project/malign)
![CI](https://github.com/tresoldi/malign/workflows/CI/badge.svg)

**MAlign** is a Python library for multiple sequence alignment with asymmetric
scoring matrices across different domains. Unlike standard alignment tools that
assume symmetric substitution costs, MAlign supports directional scoring -- the
cost of aligning symbol A with symbol B can differ from B with A.

While designed primarily for computational linguistics (e.g., historical
phonology, cognate detection), MAlign works with any hashable Python objects
and is suitable for general-purpose sequence alignment tasks.

## Key Features

- **Asymmetric scoring**: Direction-dependent alignment costs, with `from_substitution_counts()` factory for log-odds matrices from observed sound change frequencies
- **True multi-alignment**: N-dimensional alignment for up to 4 sequences (via YenKSP on N-dim graphs), with automatic progressive fallback for larger sets
- **Multiple algorithms**: Needleman-Wunsch (`anw`) and Yen's k-shortest paths (`yenksp`)
- **k-best alignments**: Return the top-k optimal alignments, not just the best one
- **Matrix learning**: Learn scoring matrices from cognate sets via EM or gradient descent
- **Feature-based scoring**: Build matrices from phonological feature distances (via [distfeat](https://github.com/tresoldi/distfeat))
- **Matrix imputation**: Fill sparse matrices using sklearn-based methods
- **Evaluation metrics**: Accuracy, precision, recall, and F1 for alignment quality

## Installation

```bash
pip install malign
```

For phonological feature-based scoring matrices:

```bash
pip install malign[features]
```

## Quick Start

### Basic Alignment

```python
import malign

alms = malign.align(["ATTCGGAT", "TACGGATTT"], k=2)
print(malign.tabulate_alms(alms))
```

### Custom Scoring Matrix

```python
matrix = malign.ScoringMatrix.from_sequences(
    sequences=[["A", "C", "G", "T"], ["A", "C", "G", "T"]],
    match=2.0, mismatch=-1.0, gap_score=-1.5,
)
alms = malign.align(["ACGT", "AGT"], k=1, matrix=matrix)
```

### Full Pipeline: Features to Evaluation

This example shows the complete workflow for linguistic alignment -- building
a scoring matrix from phonological feature distances, aligning cognate pairs,
and evaluating the results:

```python
import malign

# Build a scoring matrix from phonological feature distances
matrix = malign.ScoringMatrix.from_distfeat(
    sequences=[["n", "o", "t", "e"], ["n", "o", "tʃ", "e"]],
    gap="-", gap_score=-1.0,
)

# Align cognate sequences
alms = malign.align(
    [["n", "o", "t", "e"], ["n", "o", "tʃ", "e"]],
    k=3, matrix=matrix, method="anw",
)
print(malign.tabulate_alms(alms[:2]))

# Evaluate against gold standard
gold = malign.Alignment(
    [("n", "o", "t", "e"), ("n", "o", "tʃ", "e")], score=0.0,
)
print(f"Accuracy: {malign.alignment_accuracy(alms[0], gold):.2%}")
print(f"F1: {malign.alignment_f1(alms[0], gold):.2%}")
```

### Matrix Learning from Cognates

```python
cognate_sets = [
    [["n", "o", "t", "e"], ["n", "o", "tʃ", "e"]],
    [["f", "a", "t", "o"], ["h", "a", "d", "o"]],
]
matrix = malign.learn_matrix(cognate_sets, method="em", max_iter=10)
```

## Algorithms

| Method | Description | Best for |
|--------|-------------|----------|
| `anw` (default) | Asymmetric Needleman-Wunsch | Pairwise alignment, small k |
| `yenksp` | Yen's k-shortest paths on alignment graph | Large k, diverse alignments |
| `dumb` | Gap-padding baseline | Testing and comparison |

## Requirements

- Python >= 3.12
- numpy, scipy, scikit-learn, tabulate, PyYAML
- Optional: [distfeat](https://github.com/tresoldi/distfeat) for feature-based scoring

## Documentation

- [User Guide](docs/USER_GUIDE.md)
- [API Reference](docs/API_REFERENCE.md)
- [Algorithm Selection Guide](docs/algorithm_selection_guide.md)
- [Tutorials](docs/)

## Community Guidelines

Contributions, bug reports, and feature requests are welcome via
[GitHub issues](https://github.com/tresoldi/malign/issues) and pull requests.
See `CONTRIBUTING.md` for guidelines.

## Author and Citation

Developed by Tiago Tresoldi (tiago.tresoldi@lingfil.uu.se).

The author has received funding from the Riksbankens Jubileumsfond
(grant agreement ID: [MXM19-1087:1](https://www.rj.se/en/anslag/2019/cultural-evolution-of-texts/),
[Cultural Evolution of Texts](https://github.com/evotext/)).

During the first stages of development, the author received funding from the
European Research Council (ERC) under the European Union's Horizon 2020
research and innovation programme (grant agreement
No. [ERC Grant #715618](https://cordis.europa.eu/project/rcn/206320/factsheet/en),
[Computer-Assisted Language Comparison](https://digling.org/calc/)).

If you use `malign`, please cite it as:

  > Tresoldi, Tiago (2025). MALIGN, a library for multiple asymmetric alignments on
  > different domains. Version 0.5. Uppsala: Uppsala Universitet.

In BibTeX:

```bibtex
@misc{Tresoldi2025malign,
  author = {Tresoldi, Tiago},
  title = {MALIGN, a library for multiple asymmetric alignments on different domains. Version 0.5},
  howpublished = {\url{https://github.com/tresoldi/malign}},
  address = {Uppsala},
  publisher = {Uppsala Universitet},
  year = {2025},
}
```

## License

MIT License. See [LICENSE](LICENSE) for details.
