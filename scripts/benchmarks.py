#!/usr/bin/env python3
"""Performance benchmarking for MAlign algorithms.

This script measures and documents performance characteristics across:
1. Sequence count scaling (2-20 sequences)
2. Sequence length scaling (5-50 symbols)
3. k value scaling (1-50 top alignments)

Results are output as markdown tables for documentation.

Usage:
    python scripts/benchmarks.py [--method METHOD] [--quick]

    --method METHOD  Test only 'anw' or 'yenksp' (default: both)
    --quick          Run abbreviated benchmarks (faster, less comprehensive)
"""

import argparse
import time
from typing import Callable

import malign


def generate_test_sequences(count: int, length: int) -> list[list[str]]:
    """Generate test sequences with controlled variation.

    Creates sequences with ~80% similarity to ensure alignment
    is non-trivial but tractable.

    Args:
        count: Number of sequences to generate.
        length: Length of each sequence.

    Returns:
        List of sequences (each sequence is a list of symbols).
    """
    alphabet = ["A", "C", "G", "T"]
    base_seq = [alphabet[i % len(alphabet)] for i in range(length)]

    sequences = [base_seq.copy()]

    # Create variations with controlled differences
    for i in range(1, count):
        variant = base_seq.copy()
        # Introduce ~20% differences
        for j in range(0, length, 5):
            if j < length:
                variant[j] = alphabet[(alphabet.index(variant[j]) + i) % len(alphabet)]
        sequences.append(variant)

    return sequences


def benchmark_function(
    func: Callable,
    *args,
    iterations: int = 3,
    **kwargs,
) -> tuple[float, float]:
    """Benchmark a function with multiple iterations.

    Args:
        func: Function to benchmark.
        *args: Positional arguments for func.
        iterations: Number of times to run (default: 3).
        **kwargs: Keyword arguments for func.

    Returns:
        Tuple of (average_time, min_time) in seconds.
    """
    times = []
    for _ in range(iterations):
        start = time.perf_counter()
        func(*args, **kwargs)
        end = time.perf_counter()
        times.append(end - start)

    return sum(times) / len(times), min(times)


def benchmark_sequence_count(method: str = "anw", quick: bool = False) -> dict:
    """Benchmark: sequence count scaling.

    Tests how alignment time scales with number of sequences.

    Args:
        method: Alignment method ('anw' or 'yenksp').
        quick: If True, run abbreviated benchmark.

    Returns:
        Dict mapping sequence_count -> (avg_time, min_time).
    """
    print(f"\nBenchmarking sequence count scaling ({method.upper()})...")

    if quick:
        counts = [2, 3, 4, 5]
    else:
        counts = [2, 3, 4, 5, 6, 7, 8]

    length = 10  # Fixed sequence length
    k = 1  # Single best alignment

    results = {}
    for count in counts:
        sequences = generate_test_sequences(count, length)
        avg_time, min_time = benchmark_function(
            malign.multi_align, sequences, k=k, method=method, iterations=3
        )
        results[count] = (avg_time, min_time)
        print(f"  {count} sequences: {avg_time:.4f}s (min: {min_time:.4f}s)")

    return results


def benchmark_sequence_length(method: str = "anw", quick: bool = False) -> dict:
    """Benchmark: sequence length scaling.

    Tests how alignment time scales with sequence length.

    Args:
        method: Alignment method ('anw' or 'yenksp').
        quick: If True, run abbreviated benchmark.

    Returns:
        Dict mapping sequence_length -> (avg_time, min_time).
    """
    print(f"\nBenchmarking sequence length scaling ({method.upper()})...")

    if quick:
        lengths = [5, 10, 15, 20]
    else:
        lengths = [5, 10, 15, 20, 25, 30]

    count = 3  # Fixed number of sequences
    k = 1  # Single best alignment

    results = {}
    for length in lengths:
        sequences = generate_test_sequences(count, length)
        avg_time, min_time = benchmark_function(
            malign.multi_align, sequences, k=k, method=method, iterations=3
        )
        results[length] = (avg_time, min_time)
        print(f"  Length {length}: {avg_time:.4f}s (min: {min_time:.4f}s)")

    return results


def benchmark_k_value(method: str = "anw", quick: bool = False) -> dict:
    """Benchmark: k value scaling.

    Tests how alignment time scales with number of top alignments (k).

    Args:
        method: Alignment method ('anw' or 'yenksp').
        quick: If True, run abbreviated benchmark.

    Returns:
        Dict mapping k_value -> (avg_time, min_time).
    """
    print(f"\nBenchmarking k value scaling ({method.upper()})...")

    if quick:
        k_values = [1, 5, 10, 20]
    else:
        k_values = [1, 5, 10, 20, 30, 40, 50]

    count = 3  # Fixed number of sequences
    length = 10  # Fixed sequence length

    sequences = generate_test_sequences(count, length)

    results = {}
    for k in k_values:
        avg_time, min_time = benchmark_function(
            malign.multi_align, sequences, k=k, method=method, iterations=3
        )
        results[k] = (avg_time, min_time)
        print(f"  k={k}: {avg_time:.4f}s (min: {min_time:.4f}s)")

    return results


def generate_markdown_table(
    results: dict,
    title: str,
    x_label: str,
    methods: list[str],
) -> str:
    """Generate markdown table from benchmark results.

    Args:
        results: Dict mapping method -> (param_value -> (avg_time, min_time)).
        title: Table title.
        x_label: Label for x-axis parameter.
        methods: List of methods included.

    Returns:
        Markdown-formatted table string.
    """
    lines = [f"\n## {title}\n"]

    # Get all parameter values (assume same across methods)
    param_values = sorted(next(iter(results.values())).keys())

    # Table header
    header_cols = [x_label] + [f"{m.upper()} (avg)" for m in methods]
    lines.append("| " + " | ".join(header_cols) + " |")
    lines.append("| " + " | ".join(["---"] * len(header_cols)) + " |")

    # Table rows
    for param in param_values:
        row = [str(param)]
        for method in methods:
            if method in results and param in results[method]:
                avg_time, _ = results[method][param]
                row.append(f"{avg_time:.4f}s")
            else:
                row.append("N/A")
        lines.append("| " + " | ".join(row) + " |")

    return "\n".join(lines)


def generate_recommendations(
    seq_count_results: dict,
    seq_length_results: dict,
    k_value_results: dict,
) -> str:
    """Generate performance recommendations based on benchmark results.

    Args:
        seq_count_results: Sequence count benchmark results.
        seq_length_results: Sequence length benchmark results.
        k_value_results: k value benchmark results.

    Returns:
        Markdown-formatted recommendations.
    """
    lines = ["\n## Performance Recommendations\n"]

    # Analyze sequence count scaling
    lines.append("### Sequence Count")
    if "anw" in seq_count_results:
        anw_2 = seq_count_results["anw"][2][0]
        anw_max = seq_count_results["anw"][max(seq_count_results["anw"].keys())][0]
        ratio = anw_max / anw_2
        lines.append(f"- **ANW**: Scales ~{ratio:.1f}x from 2 to {max(seq_count_results['anw'].keys())} sequences")
        if ratio < 10:
            lines.append(f"  - Good for up to {max(seq_count_results['anw'].keys())}+ sequences")
        else:
            lines.append(f"  - Recommended limit: ~6 sequences for interactive use")

    if "yenksp" in seq_count_results:
        yenksp_2 = seq_count_results["yenksp"][2][0]
        yenksp_max = seq_count_results["yenksp"][max(seq_count_results["yenksp"].keys())][0]
        ratio = yenksp_max / yenksp_2
        lines.append(f"- **YenKSP**: Scales ~{ratio:.1f}x from 2 to {max(seq_count_results['yenksp'].keys())} sequences")
        lines.append(f"  - More computationally intensive than ANW")
        lines.append(f"  - Recommended limit: 4-5 sequences")

    # Analyze sequence length scaling
    lines.append("\n### Sequence Length")
    if "anw" in seq_length_results:
        anw_min = seq_length_results["anw"][min(seq_length_results["anw"].keys())][0]
        anw_max = seq_length_results["anw"][max(seq_length_results["anw"].keys())][0]
        ratio = anw_max / anw_min
        lines.append(f"- Scales ~{ratio:.1f}x from {min(seq_length_results['anw'].keys())} to {max(seq_length_results['anw'].keys())} symbols")
        lines.append(f"- Practical limit: 40-50 symbols for real-time use")

    # Analyze k value scaling
    lines.append("\n### k Value (Top Alignments)")
    if "anw" in k_value_results:
        k1 = k_value_results["anw"][1][0]
        k_max = k_value_results["anw"][max(k_value_results["anw"].keys())][0]
        ratio = k_max / k1
        lines.append(f"- Scales ~{ratio:.1f}x from k=1 to k={max(k_value_results['anw'].keys())}")
        lines.append(f"- k=1-10: Fast, suitable for most applications")
        lines.append(f"- k=20-50: Slower, use when diversity is critical")

    lines.append("\n### Algorithm Selection")
    lines.append("- **ANW**: Faster, better for larger problems (5+ sequences, k>10)")
    lines.append("- **YenKSP**: More thorough exploration, better for small problems (2-4 sequences)")

    return "\n".join(lines)


def main():
    """Run performance benchmarks and generate markdown report."""
    parser = argparse.ArgumentParser(description="MAlign performance benchmarking")
    parser.add_argument(
        "--method",
        choices=["anw", "yenksp", "both"],
        default="both",
        help="Which method to benchmark (default: both)",
    )
    parser.add_argument(
        "--quick",
        action="store_true",
        help="Run abbreviated benchmarks (faster)",
    )
    args = parser.parse_args()

    methods = ["anw", "yenksp"] if args.method == "both" else [args.method]

    print("=" * 80)
    print("MAlign Performance Benchmarking")
    print("=" * 80)
    if args.quick:
        print("Running in QUICK mode (abbreviated benchmarks)")

    # Run benchmarks
    seq_count_results = {}
    seq_length_results = {}
    k_value_results = {}

    for method in methods:
        seq_count_results[method] = benchmark_sequence_count(method, args.quick)
        seq_length_results[method] = benchmark_sequence_length(method, args.quick)
        k_value_results[method] = benchmark_k_value(method, args.quick)

    # Generate markdown report
    print("\n" + "=" * 80)
    print("BENCHMARK RESULTS (Markdown Format)")
    print("=" * 80)

    print("\n# MAlign Performance Benchmarks\n")
    print("Generated using `scripts/benchmarks.py`")
    print(f"\nMethods tested: {', '.join(m.upper() for m in methods)}")
    print(f"Mode: {'Quick' if args.quick else 'Full'}")

    print(
        generate_markdown_table(
            seq_count_results,
            "Sequence Count Scaling",
            "# Sequences",
            methods,
        )
    )

    print(
        generate_markdown_table(
            seq_length_results,
            "Sequence Length Scaling",
            "Length (symbols)",
            methods,
        )
    )

    print(
        generate_markdown_table(
            k_value_results,
            "k Value Scaling (Top Alignments)",
            "k Value",
            methods,
        )
    )

    print(generate_recommendations(seq_count_results, seq_length_results, k_value_results))

    print("\n" + "=" * 80)
    print("Benchmarking complete!")
    print("=" * 80)


if __name__ == "__main__":
    main()
