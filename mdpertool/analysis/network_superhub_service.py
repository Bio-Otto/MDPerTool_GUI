"""Service for identifying critical hub residues in response networks."""

from __future__ import annotations

from typing import List, Tuple

import networkx as nx


def build_superhub_rows(
    graph: nx.Graph,
    top_k: int = 5,
) -> List[Tuple[str, int, float, float, str]]:
    """
    Identify top K hub residues by betweenness centrality in the given graph.

    Parameters
    ----------
    graph : nx.Graph
        Network to analyze (typically core/intersection network)
    top_k : int
        Number of top hubs to return (default: 5)

    Returns
    -------
    List[Tuple[str, int, float, float, str]]
        Rows for superhub table.
        Each row: (node_label, degree, betweenness_centrality, closeness_centrality, connected_nodes_list)
    """
    if graph.number_of_nodes() == 0:
        return []

    # Compute centrality metrics
    betweenness = nx.betweenness_centrality(graph)
    degree_dict = dict(graph.degree())
    closeness = nx.closeness_centrality(graph)

    # Sort by betweenness centrality descending
    ranked_nodes = sorted(
        betweenness.items(),
        key=lambda x: x[1],
        reverse=True,
    )

    rows = []
    for node_label, bc_value in ranked_nodes[:top_k]:
        node_degree = degree_dict.get(node_label, 0)
        node_closeness = closeness.get(node_label, 0.0)
        
        # Get connected neighbors
        neighbors = list(graph.neighbors(node_label))
        neighbors_str = ", ".join(neighbors) if neighbors else "None"

        rows.append(
            (
                str(node_label),
                int(node_degree),
                float(bc_value),
                float(node_closeness),
                neighbors_str,
            )
        )

    return rows
