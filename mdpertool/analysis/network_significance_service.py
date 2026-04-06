"""Lightweight significance helpers for network core enrichment."""

from __future__ import annotations

from math import comb
from typing import List, Sequence, Tuple

import networkx as nx


SignificanceRow = Tuple[str, int, int, int, float, str, str, str]


def _hypergeometric_tail(population_size: int, success_population: int, sample_size: int, observed_successes: int) -> float:
    """Return P(X >= observed_successes) for the hypergeometric distribution."""
    if population_size <= 0 or sample_size <= 0 or success_population <= 0:
        return 1.0

    upper_bound = min(success_population, sample_size)
    if observed_successes > upper_bound:
        return 0.0

    denominator = comb(population_size, sample_size)
    if denominator == 0:
        return 1.0

    tail_probability = 0.0
    for k in range(observed_successes, upper_bound + 1):
        tail_probability += comb(success_population, k) * comb(population_size - success_population, sample_size - k) / denominator

    return min(1.0, max(0.0, tail_probability))


def build_network_significance_rows(
    initial_graph: nx.Graph | nx.DiGraph,
    core_graph: nx.Graph | nx.DiGraph,
    betweenness_rank_size: int = 20,
) -> List[SignificanceRow]:
    """Build a compact significance summary for the core network.

    The row reports whether the core/intersection network is enriched for the
    top-betweenness residues from the initial network.
    """
    if initial_graph is None or core_graph is None:
        return []

    initial_nodes = list(initial_graph.nodes())
    core_nodes = set(core_graph.nodes())
    if not initial_nodes or not core_nodes:
        return []

    if isinstance(initial_graph, nx.DiGraph):
        betweenness_graph = initial_graph.to_undirected()
    else:
        betweenness_graph = initial_graph

    betweenness = nx.betweenness_centrality(betweenness_graph)
    ranked_nodes = sorted(betweenness.items(), key=lambda item: item[1], reverse=True)
    top_rank_count = min(max(1, betweenness_rank_size), len(ranked_nodes))
    top_nodes = [node for node, _ in ranked_nodes[:top_rank_count]]

    overlap_nodes = [node for node in top_nodes if node in core_nodes]
    overlap_count = len(overlap_nodes)
    top_nodes_text = ", ".join(top_nodes) if top_nodes else "None"
    overlap_nodes_text = ", ".join(overlap_nodes) if overlap_nodes else "None"

    p_value = _hypergeometric_tail(
        population_size=len(initial_nodes),
        success_population=len(top_nodes),
        sample_size=len(core_nodes),
        observed_successes=overlap_count,
    )

    return [
        (
            f"Core enrichment (Top-{top_rank_count} candidates)",
            len(core_nodes),
            len(top_nodes),
            overlap_count,
            p_value,
            top_nodes_text,
            overlap_nodes_text,
            f"Compared {top_rank_count} top-betweenness candidates against the core network; overlap = {overlap_count}",
        )
    ]
