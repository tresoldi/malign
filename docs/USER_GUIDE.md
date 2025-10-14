# MAlign User Guide

**Version**: 0.4.0 (in development)

**Note**: Comprehensive user guide to be completed in Phase 4. This is a structural template.

---

## Table of Contents

1. [Getting Started](#getting-started)
2. [Core Concepts](#core-concepts)
3. [Matrix Management](#matrix-management)
4. [Algorithms](#algorithms)
5. [Advanced Topics](#advanced-topics)
6. [Performance & Limits](#performance--limits)
7. [Troubleshooting](#troubleshooting)
8. [API Reference](#api-reference)

---

## 1. Getting Started

### Installation

```bash
pip install malign
```

### Quick Start

```python
import malign

# Basic alignment
sequences = ["ATTCG", "ATGC"]
alignments = malign.align(sequences, k=2)

# View results
print(malign.tabulate_alms(alignments))
```

**TODO**: Expand with more examples in Phase 4.

---

## 2. Core Concepts

### Asymmetric Scoring

**TODO**: Explain why asymmetry matters for linguistics in Phase 4.

### True Multi-Alignment

**TODO**: Compare with progressive methods (UPGMA, lingpy) in Phase 4.

### K-Best Results

**TODO**: Explain interpretation and use cases in Phase 4.

---

## 3. Matrix Management

### Construction

**TODO**: Document all builder methods in Phase 4.

### Serialization (YAML)

**TODO**: Document YAML format with examples in Phase 4.

### Optimization

**TODO**: Document matrix learning in Phase 4.

---

## 4. Algorithms

### ANW (Asymmetric Needleman-Wunsch)

- **Best for**: TODO
- **Complexity**: TODO
- **When to use**: TODO

### YenKSP (Yen's K-Shortest Paths)

- **Best for**: TODO
- **Complexity**: TODO
- **When to use**: TODO

**TODO**: Add decision flowchart in Phase 4 based on Phase 3 benchmarks.

---

## 5. Advanced Topics

### Batch Processing

**TODO**: Document align_batch() in Phase 4.

### Validation Metrics

**TODO**: Document all metrics in Phase 4.

### Integration with External Libraries

- **freqprob**: TODO
- **asymcat**: TODO
- **nhandu**: TODO

---

## 6. Performance & Limits

### Computational Limits

**TODO**: Document based on Phase 3 benchmarks.

- Recommended sequence count: TBD
- Maximum practical: TBD
- Performance tips: TBD

### Benchmarks

**TODO**: Include ANW vs YenKSP comparison from Phase 3.

---

## 7. Troubleshooting

### Common Issues

**TODO**: Populate based on user feedback.

### Performance Issues

**TODO**: Add profiling and optimization tips.

---

## 8. API Reference

See [API_REFERENCE.md](API_REFERENCE.md) for complete API documentation.

---

**Completion Status**: Phase 1 template created. Content to be added in Phase 4.
