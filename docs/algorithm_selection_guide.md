# Algorithm Selection Guide

This guide helps you choose between alignment methods and understand how MAlign dispatches alignment strategies automatically.

## Automatic Multi-Sequence Dispatch

MAlign automatically selects the best alignment strategy based on input size:

```
┌─────────────────────────────────────┐
│  How many sequences?                │
└──────────────┬──────────────────────┘
               │
      ┌────────┴────────┐
      │                 │
    N = 2             N >= 3
      │                 │
      v                 v
┌───────────┐    ┌──────────────────┐
│  Direct   │    │ Grid small       │
│  pairwise │    │ enough for       │
│  (NW or   │    │ N-dim? (N<=4,    │
│  YenKSP)  │    │ cells<=200K)     │
└───────────┘    └────────┬─────────┘
                    ┌─────┴─────┐
                   Yes         No
                    │           │
                    v           v
              ┌──────────┐ ┌──────────────┐
              │ True     │ │ UPGMA        │
              │ N-dim    │ │ progressive  │
              │ YenKSP   │ │ alignment    │
              │ (optimal)│ │ (beam search │
              └──────────┘ │ for k-best)  │
                           └──────────────┘
```

This dispatch is transparent -- just call `malign.align()` and the right strategy is selected.

**N-dimensional alignment** (N=3-4, small grids) finds the globally optimal alignment across all sequences simultaneously, not a pairwise approximation.

**UPGMA progressive alignment** (N=5+ or large grids) builds a guide tree from pairwise distances, then merges sequences along the tree with beam search for k-best results.

## ANW vs YenKSP

Within each dispatch tier, you choose between two alignment algorithms:

### Decision Flowchart

```
┌─────────────────────────────────────┐
│  How many sequences do you have?    │
└──────────────┬──────────────────────┘
               │
      ┌────────┴────────┐
      │                 │
    <= 4               > 4
      │                 │
      v                 v
┌─────────────┐    ┌──────────┐
│ Do you need │    │ Use ANW  │
│ k > 10?     │    │          │
└─────┬───────┘    │ (YenKSP  │
      │            │  too     │
   ┌──┴──┐         │  slow)   │
   │     │         └──────────┘
  Yes   No
   │     │
   v     v
┌─────┐ ┌─────────────┐
│ ANW │ │   YenKSP    │
│     │ │             │
│ (k  │ │ (thorough   │
│ >10)│ │  for small  │
└─────┘ │  problems)  │
        └─────────────┘
```

**Quick Decision Rules:**
- **2-4 sequences, k <= 10**: Use **YenKSP** (more thorough)
- **2-4 sequences, k > 10**: Use **ANW** (faster for large k)
- **5+ sequences**: Use **ANW** (YenKSP too slow)
- **Real-time/interactive**: Use **ANW** (generally faster)
- **Maximum quality, small problem**: Use **YenKSP**

## Comparison Table

| Feature | ANW | YenKSP |
|---------|-----|--------|
| **Speed** | Faster (especially for k > 10) | Slower |
| **Quality** | High quality, near-optimal | Highest quality, exhaustive search |
| **Sequence Count** | Good for 5-8+ sequences | Best for 2-4 sequences |
| **k Value Scaling** | Near-linear (~2x for k=1->20) | More expensive for large k |
| **Use Cases** | Large problems, interactive use, k > 10 | Small problems, maximum quality |
| **Algorithm** | A* search + NW alignment | Yen's k-shortest paths |

## Detailed Pros & Cons

### ANW (A* + Needleman-Wunsch)

**Pros:**
- Faster, especially for k > 10 and many sequences
- Scales better: can handle 5-8 sequences reasonably
- Near-linear k scaling: 2x slower for k=20 vs k=1
- Fast enough for real-time/interactive use
- Heuristic guidance: A* search efficiently explores space

**Cons:**
- Heuristic-based: may miss some optimal solutions
- Slightly lower quality than YenKSP for small problems

**Best For:**
- Problems with 5+ sequences
- High k values (k > 10)
- Interactive/real-time applications
- When speed matters more than absolute optimality

### YenKSP (Yen's k-Shortest Paths)

**Pros:**
- Exhaustive: guarantees finding true k-shortest paths
- Highest alignment quality
- Deterministic: always finds optimal k-best alignments
- Strong theoretical guarantees

**Cons:**
- Slower, especially for large k and many sequences
- 362x slower for 5 vs 2 sequences
- Impractical for 5+ sequences
- Becomes very slow for k > 20

**Best For:**
- Problems with 2-4 sequences
- Small to moderate k (k <= 10)
- When maximum quality is critical
- Batch processing (not interactive)
- Research/benchmarking

## Post-Processing: Block Merging

After alignment (with either method), you can apply block detection
to merge complementary-gap columns (diphthongization, metathesis):

```python
# Automatic: merge during alignment
alms = malign.align(sequences, k=1, merge_blocks=True, max_block_size=2)

# Manual: merge an existing alignment
merged = malign.merge_alignment_blocks(alignment, max_block_size=2)
```

Block merging is post-processing only -- it does not affect the alignment
scoring or algorithm choice.

## Performance Benchmarks

Based on `scripts/benchmarks.py` results:

### Sequence Count Impact (k=1)

| Sequences | ANW Time | Notes |
|-----------|----------|-------|
| 2 | 0.001s | Instant |
| 3 | 0.004s | Very fast |
| 4 | 0.053s | Fast |
| 5 | 0.397s | Usable |
| 6-8 | 1-10s (est) | Batch only |

**Scaling**: ~362x from 2->5 sequences (exponential)

### k Value Impact (3 sequences)

| k | ANW Time | Use Case |
|---|----------|----------|
| 1 | 0.004s | Single best |
| 5 | 0.004s | Top 5 |
| 10 | 0.006s | Top 10 |
| 20 | 0.008s | Diversity |

**Scaling**: ~2x from k=1->20 (nearly linear)

### Sequence Length Impact (3 sequences, k=1)

| Length | Time | Notes |
|--------|------|-------|
| 5 | 0.004s | Short |
| 10 | 0.004s | Medium |
| 15 | 0.005s | Long |
| 20 | 0.010s | Very long |

**Scaling**: ~2.7x from 5->20 symbols (sub-quadratic)

## Practical Recommendations

### Interactive Applications
Use **ANW** with conservative parameters:
- k <= 10 for instant results
- <= 5 sequences for real-time
- <= 20 symbols per sequence

### Batch Processing
Either method works:
- **ANW**: For large problems (5+ sequences, k > 10)
- **YenKSP**: For maximum quality on small problems

### Research/Benchmarking
Use **YenKSP** for ground truth:
- Guarantees true k-best alignments
- Use as gold standard for evaluating other methods
- Limited to small problems (<=4 sequences)

### Production Systems
Use **ANW** for reliability:
- Predictable performance
- Handles varied input sizes
- Good quality/speed trade-off

## Examples

### Example 1: Aligning 3 cognate words
```python
# Small problem, want top 10 alignments
sequences = [["k", "a", "t"], ["c", "a", "t"], ["k", "a", "t", "z"]]

# Use YenKSP for maximum quality
alms = malign.align(sequences, k=10, method="yenksp")
```

### Example 2: Aligning 6 language forms
```python
# Larger problem -- automatically uses progressive alignment
sequences = [[...], [...], [...], [...], [...], [...]]  # 6 sequences

# Use ANW (YenKSP would be too slow)
alms = malign.align(sequences, k=5, method="anw")
```

### Example 3: Exploring alignment space
```python
# Want to see many alternatives (k=50)
sequences = [["A", "B", "C"], ["A", "B", "D"]]

# Use ANW (k=50 would be slow in YenKSP)
alms = malign.align(sequences, k=50, method="anw")
```

### Example 4: Block-aware alignment
```python
# Diphthongization pattern
alms = malign.align(
    [["a"], ["j", "e"]],
    k=1,
    merge_blocks=True,
)
# Sequence 2 gets compound symbol ("j", "e")
```

## Summary

**Default Choice**: Use **ANW** unless you have a specific reason to use YenKSP.

**Use YenKSP** when:
- You have <= 4 sequences
- You need k <= 10 alignments
- Maximum quality is critical
- You're doing research/benchmarking

**Use ANW** when:
- You have 5+ sequences
- You need k > 10 alignments
- Speed matters (interactive use)
- You're building a production system
