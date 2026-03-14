# malign examples

Runnable scripts showcasing the `malign` API on real linguistic and biological datasets. Each script loads data, runs alignment, writes output files to `output/`, and prints a summary.

## Scripts

| # | Script | Description | API features | Dataset | Extra deps |
|---|--------|-------------|-------------|---------|------------|
| 01 | `01_dna_alignment.py` | Basic alignment, k-best, method comparison | `align()`, `tabulate_alms()`, `DNA_MATRIX` | Hardcoded DNA sequences | -- |
| 02 | `02_bdpa_evaluation.py` | Gold-standard evaluation of predicted alignments | `align()`, `Alignment()`, `alignment_accuracy()`, `alignment_precision_recall()`, `alignment_f1()` | `data/bdpa/romance.tsv` | -- |
| 03 | `03_northeuralex_pipeline.py` | distfeat -> bootstrap -> learn pipeline | `ScoringMatrix.from_distfeat()`, `bootstrap_matrix()`, `learn_matrix()` | `data/northeuralex/italic_by_cogid.tsv` | `pip install malign[features]` |
| 04 | `04_cmudict_cross_domain.py` | Cross-domain alignment (letters <-> phonemes) | `align()`, `ScoringMatrix.from_substitution_counts()` | `data/cmudict/cmudict.tsv` | -- |
| 05 | `05_wiktionary_multilingual.py` | Romance cognate alignment (2-way and 3-way) | `align()`, `bootstrap_matrix()` | `data/wiktionary/wiktionary.tsv` | -- |
| 06 | `06_neurodecipher_scripts.py` | Disjoint-domain bootstrap (Greek alphabet vs Linear B syllabary) | `align()`, `bootstrap_matrix()` | `data/neurodecipher/greek_linearb.tsv` | -- |
| 07 | `07_lexibank_learning.py` | Large-scale matrix learning, gold evaluation | `learn_matrix()`, `align()`, `alignment_accuracy()` | `data/lexibank/forms.tsv` | -- |

## Running

```bash
# Run a single example
python examples/01_dna_alignment.py

# Run all (script 03 requires malign[features])
for f in examples/0*.py; do python "$f"; done
```

Output files are written to `examples/output/`.

## Data

All datasets live under `examples/data/`. See each subdirectory's `README.md` for provenance and licensing.
