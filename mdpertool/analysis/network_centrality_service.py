"""Centrality helpers for responsive exact network analysis workflows."""

from __future__ import annotations

from typing import Any, Dict

import networkx as nx

from networkx.algorithms.centrality.betweenness import (  # type: ignore
    _accumulate_basic,
    _accumulate_endpoints,
    _rescale,
    _single_source_dijkstra_path_basic,
    _single_source_shortest_path_basic,
)


def calculate_betweenness_with_progress(
    graph: nx.Graph,
    normalized: bool = True,
    weight: Any = None,
    endpoints: bool = False,
    progress_callback: Any = None,
) -> Dict[Any, float]:
    """Compute exact betweenness while emitting coarse-grained progress updates."""
    if graph is None or graph.number_of_nodes() == 0:
        return {}

    betweenness: Dict[Any, float] = dict.fromkeys(graph, 0.0)
    nodes = list(graph)
    total_nodes = len(nodes)

    for index, source in enumerate(nodes, start=1):
        if weight is None:
            stack, predecessors, sigma, _ = _single_source_shortest_path_basic(graph, source)
        else:
            stack, predecessors, sigma, _ = _single_source_dijkstra_path_basic(graph, source, weight)

        if endpoints:
            betweenness, _ = _accumulate_endpoints(betweenness, stack, predecessors, sigma, source)
        else:
            betweenness, _ = _accumulate_basic(betweenness, stack, predecessors, sigma, source)

        if callable(progress_callback) and (index == 1 or index == total_nodes or index % 5 == 0):
            progress_callback(index, total_nodes)

    # NetworkX private helper signatures vary between releases.
    try:
        betweenness = _rescale(
            betweenness,
            len(graph),
            normalized=normalized,
            directed=graph.is_directed(),
            k=None,
            endpoints=endpoints,
        )
    except TypeError:
        try:
            betweenness = _rescale(
                betweenness,
                len(graph),
                normalized=normalized,
                directed=graph.is_directed(),
                endpoints=endpoints,
            )
        except TypeError:
            betweenness = _rescale(
                betweenness,
                len(graph),
                normalized,
                graph.is_directed(),
                endpoints,
            )
    return betweenness
