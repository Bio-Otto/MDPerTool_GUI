"""Motif summary helpers for directed network analysis."""

from __future__ import annotations

from itertools import combinations, permutations
from typing import Dict, List, Sequence, Tuple

import networkx as nx


MotifRow = Tuple[str, str, int, int, str, str]


def _canonical_motif_key(subgraph: nx.DiGraph) -> str:
    """Return canonical motif key by permuting adjacency matrices.

    Subgraphs are small (3 or 4 nodes), so brute-force canonicalization is practical.
    """
    nodes = list(subgraph.nodes())
    size = len(nodes)
    best_key = None

    for perm in permutations(range(size)):
        bit_list = []
        for i in perm:
            u = nodes[i]
            for j in perm:
                v = nodes[j]
                bit_list.append('1' if subgraph.has_edge(u, v) else '0')
        key = ''.join(bit_list)
        if best_key is None or key < best_key:
            best_key = key

    return best_key or ''


def _count_connected_induced_motifs(
    graph: nx.DiGraph,
    motif_size: int,
    max_combinations: int,
) -> Tuple[Dict[str, int], int, str]:
    """Count connected induced motifs for a specific size.

    Returns motif count map, total connected motif count, and status text.
    """
    node_count = graph.number_of_nodes()
    if node_count < motif_size:
        return {}, 0, "not_enough_nodes"

    combo_estimate = 1
    for i in range(motif_size):
        combo_estimate = combo_estimate * (node_count - i) // (i + 1)
    if combo_estimate > max_combinations:
        return {}, 0, f"skipped_too_many_combinations({combo_estimate})"

    motif_counter: Dict[str, int] = {}
    connected_count = 0

    for node_group in combinations(graph.nodes(), motif_size):
        subgraph = nx.DiGraph(graph.subgraph(node_group).copy())
        if subgraph.number_of_edges() == 0:
            continue
        if not nx.is_weakly_connected(subgraph):
            continue

        connected_count += 1
        key = _canonical_motif_key(subgraph)
        motif_counter[key] = motif_counter.get(key, 0) + 1

    return motif_counter, connected_count, "ok"


def build_motif_summary_rows(
    graph: nx.DiGraph,
    scope_name: str,
    motif_sizes: Sequence[int] = (3, 4),
    max_combinations: int = 400000,
    top_k_per_size: int = 10,
) -> List[MotifRow]:
    """Build motif summary rows for analysis table.

    Row schema:
    (size, motif_id, edge_count, occurrence, frequency_pct, scope)
    """
    if graph is None:
        return []

    rows: List[MotifRow] = []

    for motif_size in motif_sizes:
        motif_counter, connected_count, status = _count_connected_induced_motifs(
            graph=graph,
            motif_size=motif_size,
            max_combinations=max_combinations,
        )

        if status != "ok":
            rows.append(
                (
                    f"{motif_size}-node",
                    "N/A",
                    0,
                    0,
                    status,
                    scope_name,
                )
            )
            continue

        if connected_count == 0:
            rows.append(
                (
                    f"{motif_size}-node",
                    "N/A",
                    0,
                    0,
                    "0.00%",
                    scope_name,
                )
            )
            continue

        sorted_items = sorted(motif_counter.items(), key=lambda item: item[1], reverse=True)[:top_k_per_size]
        for idx, (motif_key, occurrence) in enumerate(sorted_items, start=1):
            edge_count = motif_key.count('1')
            frequency_pct = (occurrence / connected_count) * 100.0
            rows.append(
                (
                    f"{motif_size}-node",
                    f"M{motif_size}_{idx:02d}",
                    edge_count,
                    occurrence,
                    f"{frequency_pct:.2f}%",
                    scope_name,
                )
            )

    return rows
