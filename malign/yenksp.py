"""Module with functions for the k-best pairwise alignment.

Uses a pure Python implementation of Yen's K-Shortest Paths algorithm
with a simple adjacency-dict graph representation.
"""

import heapq
from collections.abc import Hashable, Sequence
from itertools import pairwise

from .alignment import Alignment
from .scoring_matrix import ScoringMatrix
from .utils import pad_sequence, score_alignment, sort_alignments

# Type alias for graph: node -> {neighbor -> weight}
type Node = tuple[int, int]
type Graph = dict[Node, dict[Node, float]]


def _new_graph() -> Graph:
    """Create an empty graph."""
    return {}


def _add_edge(graph: Graph, u: Node, v: Node, weight: float) -> None:
    """Add a directed edge from u to v with the given weight."""
    if u not in graph:
        graph[u] = {}
    graph[u][v] = weight
    # Ensure v exists as a node even if it has no outgoing edges
    if v not in graph:
        graph[v] = {}


def _neighbors(graph: Graph, node: Node) -> dict[Node, float]:
    """Return neighbors and weights for a node."""
    return graph.get(node, {})


def _dijkstra(
    graph: Graph,
    source: Node,
    target: Node,
    excluded_nodes: set[Node] | None = None,
    excluded_edges: set[tuple[Node, Node]] | None = None,
) -> list[Node] | None:
    """Find shortest path from source to target using Dijkstra's algorithm.

    Args:
        graph: Adjacency-dict graph.
        source: Start node.
        target: End node.
        excluded_nodes: Nodes to ignore during search.
        excluded_edges: Edges to ignore during search.

    Returns:
        Shortest path as list of nodes, or None if no path exists.
    """
    if excluded_nodes and source in excluded_nodes:
        return None

    dist: dict[Node, float] = {source: 0.0}
    prev: dict[Node, Node] = {}
    heap: list[tuple[float, int, Node]] = [(0.0, 0, source)]
    counter = 1  # Tie-breaker for heap

    while heap:
        d, _, u = heapq.heappop(heap)

        if u == target:
            # Reconstruct path
            path = [target]
            while path[-1] != source:
                path.append(prev[path[-1]])
            path.reverse()
            return path

        if d > dist.get(u, float("inf")):
            continue

        for v, w in _neighbors(graph, u).items():
            if excluded_nodes and v in excluded_nodes:
                continue
            if excluded_edges and (u, v) in excluded_edges:
                continue

            new_dist = d + w
            if new_dist < dist.get(v, float("inf")):
                dist[v] = new_dist
                prev[v] = u
                heapq.heappush(heap, (new_dist, counter, v))
                counter += 1

    return None


def _path_cost(graph: Graph, path: list[Node]) -> float:
    """Compute total edge-weight cost of a path."""
    cost = 0.0
    for i in range(len(path) - 1):
        cost += graph[path[i]][path[i + 1]]
    return cost


def yen_ksp(
    graph: Graph,
    source: Node,
    target: Node,
    k: int,
) -> list[list[Node]]:
    """Find k shortest simple paths using Yen's algorithm.

    Args:
        graph: Adjacency-dict graph.
        source: Start node.
        target: End node.
        k: Number of shortest paths to find.

    Returns:
        List of up to k shortest paths (each a list of nodes).
    """
    # Find the first shortest path
    first_path = _dijkstra(graph, source, target)
    if first_path is None:
        return []

    a: list[list[Node]] = [first_path]  # Confirmed shortest paths
    b: list[tuple[float, int, list[Node]]] = []  # Candidate paths (heap)
    b_counter = 0  # Tie-breaker
    b_set: set[tuple[Node, ...]] = set()  # Dedup candidates

    for i in range(1, k):
        prev_path = a[i - 1]

        for j in range(len(prev_path) - 1):
            spur_node = prev_path[j]
            root_path = prev_path[: j + 1]

            # Exclude edges that share the same root path prefix
            excluded_edges: set[tuple[Node, Node]] = set()
            for path in a:
                if len(path) > j and path[: j + 1] == root_path:
                    excluded_edges.add((path[j], path[j + 1]))

            # Exclude root path nodes (except spur node) to ensure simple paths
            excluded_nodes = set(root_path[:-1])

            # Find spur path from spur_node to target
            spur_path = _dijkstra(graph, spur_node, target, excluded_nodes, excluded_edges)

            if spur_path is not None:
                total_path = root_path[:-1] + spur_path
                path_key = tuple(total_path)
                if path_key not in b_set:
                    b_set.add(path_key)
                    cost = _path_cost(graph, total_path)
                    heapq.heappush(b, (cost, b_counter, total_path))
                    b_counter += 1

        if not b:
            break

        # Pop the best candidate
        _, _, best_path = heapq.heappop(b)
        a.append(best_path)

    return a


def compute_graph(
    seq_a: Sequence[Hashable],
    seq_b: Sequence[Hashable],
    matrix: ScoringMatrix,
) -> Graph:
    """Compute a weighted directed graph for alignment.

    Builds a weighted directed graph for a pairwise alignment between two
    sequences, analogous to a Needleman-Wunsch alignment matrix. All possible
    diagonal, vertical, and horizontal transitions are computed.

    Args:
        seq_a: First sequence to be aligned.
        seq_b: Second sequence to be aligned.
        matrix: Scoring matrix for the alignment.

    Returns:
        A directed graph as adjacency dict.
    """
    max_score = max(matrix.scores.values())

    # Add gaps to the beginning of both sequences
    seq_a = pad_sequence(seq_a, matrix.gap)
    seq_b = pad_sequence(seq_b, matrix.gap)

    graph = _new_graph()
    for i in range(len(seq_a) - 1, -1, -1):
        symbol_a = seq_a[i]

        for j in range(len(seq_b) - 1, -1, -1):
            symbol_b = seq_b[j]

            if i == 0 and j == 0:
                dig_score = None
                hor_score = None
                ver_score = None
            elif i == 0:
                dig_score = None
                hor_score = None
                # (i, j-1) -> (i, j): consume seq_b symbol, gap in seq_a
                ver_score = matrix[matrix.gap, symbol_b]
            elif j == 0:
                dig_score = None
                # (i-1, j) -> (i, j): consume seq_a symbol, gap in seq_b
                hor_score = matrix[symbol_a, matrix.gap]
                ver_score = None
            else:
                dig_score = matrix[symbol_a, symbol_b]
                hor_score = matrix[symbol_a, matrix.gap]
                ver_score = matrix[matrix.gap, symbol_b]

            if dig_score is not None:
                _add_edge(graph, (i - 1, j - 1), (i, j), max_score - dig_score)
            if hor_score is not None:
                _add_edge(graph, (i - 1, j), (i, j), max_score - hor_score)
            if ver_score is not None:
                _add_edge(graph, (i, j - 1), (i, j), max_score - ver_score)

    return graph


def build_align(
    path: list[tuple[int, int]],
    seq_a: Sequence[Hashable],
    seq_b: Sequence[Hashable],
    gap: Hashable = "-",
) -> tuple[Sequence[Hashable], Sequence[Hashable]]:
    """Build a pairwise alignment from a path of sequence indexes.

    Args:
        path: A list of (row, col) index tuples from the alignment graph.
        seq_a: First sequence (vertical border).
        seq_b: Second sequence (horizontal border).
        gap: Gap symbol. Defaults to "-".

    Returns:
        The aligned sequences for A and B.
    """
    seq_a = list(seq_a[path[0][0] :])
    seq_b = list(seq_b[path[0][1] :])

    alm_a = []
    alm_b = []
    for source, target in pairwise(path):
        if target[1] == source[1]:
            # vertical movement
            alm_a.append(seq_a.pop(0))
            alm_b.append(gap)
        elif target[0] == source[0]:
            # horizontal movement
            alm_a.append(gap)
            alm_b.append(seq_b.pop(0))
        else:
            # diagonal movement
            alm_a.append(seq_a.pop(0))
            alm_b.append(seq_b.pop(0))

    return alm_a, alm_b


def align(
    graph: Graph,
    ne_loc: tuple[int, int],
    sw_loc: tuple[int, int],
    seq_a: Sequence[Hashable],
    seq_b: Sequence[Hashable],
    matrix: ScoringMatrix,
    n_paths: int | None = None,
) -> list[Alignment]:
    """Return the k best alignments in terms of costs.

    Args:
        graph: A directed weighted graph for the alignment.
        ne_loc: Source node (top-left corner).
        sw_loc: Target node (bottom-right corner).
        seq_a: First sequence (vertical border).
        seq_b: Second sequence (horizontal border).
        matrix: Scoring matrix.
        n_paths: Number of alignment paths to collect.

    Returns:
        A sorted list of best alignments.
    """
    n_paths = n_paths or min(len(seq_a), len(seq_b))

    paths = yen_ksp(graph, ne_loc, sw_loc, n_paths)

    alignments = []
    for path in paths:
        alm_seq_a, alm_seq_b = build_align(path, seq_a, seq_b, gap=matrix.gap)
        score = score_alignment([alm_seq_a, alm_seq_b], matrix)
        alignments.append(Alignment([alm_seq_a, alm_seq_b], score))

    return sort_alignments(alignments)


def yenksp_align(
    seq_a: Sequence[Hashable],
    seq_b: Sequence[Hashable],
    k: int | None = 4,
    matrix: ScoringMatrix | None = None,
    ne_loc: tuple[int, int] = (0, 0),
    sw_loc: tuple[int, int] | None = None,
) -> list[Alignment]:
    """Perform pairwise alignment with the Yen K Shortest Paths method.

    Args:
        seq_a: The first sequence to be aligned.
        seq_b: The second sequence to be aligned.
        k: The number of paths to compute. Defaults to 4.
        matrix: The scoring matrix. If not provided, an identity matrix is created.
        ne_loc: Source node. Defaults to (0, 0).
        sw_loc: Target node. Defaults to (len(seq_a), len(seq_b)).

    Returns:
        A sorted list of best alignments.
    """
    if not matrix:
        matrix = ScoringMatrix.from_sequences([seq_a, seq_b])

    sw_loc = sw_loc or (len(seq_a), len(seq_b))

    graph = compute_graph(seq_a, seq_b, matrix)

    return align(graph, ne_loc, sw_loc, seq_a, seq_b, matrix, n_paths=k)
