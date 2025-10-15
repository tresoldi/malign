# MAlign Test Data

This directory contains gold standard test data extracted from the **Arca Verborum** project for regression testing and matrix learning evaluation.

## Directory Structure

```
tests/data/
├── README.md (this file)
└── cognates/
    ├── regression_test_set.yml
    ├── learning_training_set.yml
    ├── learning_eval_set.yml
    └── integration_examples.yml
```

## Data Source

All test data is extracted from:
- **Source**: `resources/forms_with_alignments.csv`
- **Project**: Arca Verborum - comprehensive linguistic database
- **Size**: 451,935 forms across 111,142 cognate sets from 58 datasets
- **Content**: Expert-curated phonological alignments across language families

## Dataset Descriptions

### `regression_test_set.yml`
**Purpose**: Regression testing - verify alignment quality against gold standards

- **Count**: 100 cognate sets
- **Selection Criteria**:
  - 3-10 forms per set (diverse multiwise scenarios)
  - Contains gaps (tests gap handling)
  - Diverse across linguistic families (18+ datasets)
  - High-quality gold alignments
- **Usage**: `tests/test_regression.py`
- **Metrics**: Baseline ~63% accuracy with identity matrix

**Example**:
```yaml
- id: bowernpny_28-other
  dataset: bowernpny
  gloss: OTHER
  forms:
    - language_id: bowernpny_Arabana
      segments: [ŋ, u, ɹ, u]
      alignment: [ŋ, u, ɹ, u, '-', '-', '-', '-']
```

### `learning_training_set.yml`
**Purpose**: Matrix learning training data

- **Count**: 200 cognate sets
- **Selection Criteria**:
  - Larger, more diverse sample for training
  - Same quality standards as regression set
  - No overlap with evaluation set
- **Usage**: Train EM and gradient descent algorithms
- **Expected**: Learned matrices should significantly improve alignment accuracy

### `learning_eval_set.yml`
**Purpose**: Matrix learning evaluation data

- **Count**: 100 cognate sets
- **Selection Criteria**:
  - **No overlap** with training set (clean train/test split)
  - Similar distribution to training set
  - Independent validation of learned matrices
- **Usage**: Evaluate generalization of learned matrices
- **Metrics**: Compare learned matrix accuracy vs identity matrix baseline

### `integration_examples.yml`
**Purpose**: Small diverse examples for integration tests

- **Count**: 20 cognate sets
- **Selection Criteria**:
  - Size diversity: 2-3 forms (33%), 4-6 forms (33%), 7-10 forms (33%)
  - Varied linguistic families
  - Quick smoke tests
- **Usage**: `tests/test_integration.py` for fast validation
- **Purpose**: Verify end-to-end pipelines work correctly

## Data Format

All datasets follow the same YAML schema:

```yaml
description: <dataset purpose>
source: Arca Verborum (forms_with_alignments.csv)
count: <number of cognate sets>
cognate_sets:
  - id: <unique identifier>
    dataset: <source dataset name>
    gloss: <concepticon gloss/meaning>
    parameter_id: <parameter identifier>
    forms:
      - language_id: <language identifier>
        glottolog_name: <standardized language name>
        segments: [<unaligned phonemes>]
        alignment: [<aligned phonemes with gaps>]
```

## Usage Examples

### Load Regression Data
```python
from tests.gold_data_utils import load_cognate_sets

cognate_sets = load_cognate_sets("tests/data/cognates/regression_test_set.yml")
print(f"Loaded {len(cognate_sets)} cognate sets")
```

### Convert to Sequences
```python
from tests.gold_data_utils import cognate_set_to_sequences, cognate_set_to_gold_alignment

sequences = cognate_set_to_sequences(cognate_sets[0])
gold_alignment = cognate_set_to_gold_alignment(cognate_sets[0])
```

### Filter by Criteria
```python
from tests.gold_data_utils import filter_cognate_sets

# Get only small sets (2-3 forms) with gaps
small_gappy = filter_cognate_sets(
    cognate_sets,
    min_forms=2,
    max_forms=3,
    has_gaps=True
)
```

## Test Utilities

**Module**: `tests/gold_data_utils.py`

Provides:
- `load_cognate_sets(yaml_path)`: Load YAML datasets
- `GoldCognateSet`: Named tuple for cognate set structure
- `GoldForm`: Named tuple for individual language forms
- `cognate_set_to_sequences()`: Extract unaligned sequences
- `cognate_set_to_gold_alignment()`: Extract gold alignment
- `filter_cognate_sets()`: Filter by size, dataset, gaps

## Statistics

### Dataset Coverage

| Dataset | Cognate Sets | Avg Forms/Set |
|---------|--------------|---------------|
| bowernpny | 2,430 | 4.2 |
| smithborneo | 2,141 | 3.8 |
| castroyi | 1,156 | 5.1 |
| (+ 15 more) | ... | ... |

### Form Characteristics

- **Total Forms**: 451,935 across all datasets
- **With Gaps**: 17.2% (78,000+)
- **Avg Alignment Length**: 12.3 symbols
- **Avg Sequence Length**: 8.1 symbols (before alignment)

## Data Quality

All alignments are:
- ✅ **Expert-curated**: Created by linguistic researchers
- ✅ **Multi-language**: Span diverse language families
- ✅ **Phonologically-informed**: Based on sound correspondence theory
- ✅ **Validated**: Part of published research database

## Maintenance

### Regenerating Datasets

If you need to regenerate the test datasets:

```bash
# 1. Ensure resources/forms_with_alignments.csv is present
# 2. Run extraction scripts
python scripts/analyze_gold_data.py    # Analyze source data
python scripts/extract_test_data.py    # Extract test sets
```

### Adding New Datasets

To add new test datasets:
1. Follow the YAML schema above
2. Place in `tests/data/cognates/`
3. Update this README
4. Add tests in appropriate test file
5. Update `tests/gold_data_utils.py` if new utilities needed

## References

- **Arca Verborum**: https://github.com/lexibank/arca
- **Cognate Set Alignment**: Expert manual alignments from comparative linguistics
- **Gold Standard**: Used for validation, not training (except learning sets)

## License

Test data is derived from Arca Verborum (CC-BY-4.0 licensed research data).
