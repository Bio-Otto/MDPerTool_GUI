"""Network summary helpers for Initial/Union/Intersection reporting."""

from __future__ import annotations

from typing import Iterable, List, Sequence, Tuple

import networkx as nx


SummaryRow = Tuple[str, int, int, int, int, float, str, float, float]


def _as_digraph(graph: nx.Graph | nx.DiGraph) -> nx.DiGraph:
    if isinstance(graph, nx.DiGraph):
        return graph.copy()
    return nx.DiGraph(graph)


def _largest_undirected_component(graph: nx.DiGraph) -> nx.Graph:
    if graph.number_of_nodes() == 0:
        return nx.Graph()

    undirected = graph.to_undirected()
    if undirected.number_of_nodes() == 0:
        return nx.Graph()

    components = list(nx.connected_components(undirected))
    if not components:
        return nx.Graph()
    largest_nodes = max(components, key=len)
    return undirected.subgraph(largest_nodes).copy()


def build_union_graph(graphs: Sequence[nx.Graph | nx.DiGraph]) -> nx.DiGraph:
    if not graphs:
        return nx.DiGraph()
    directed_graphs = [_as_digraph(graph) for graph in graphs]
    return nx.compose_all(directed_graphs)


def build_intersection_graph(graphs: Sequence[nx.Graph | nx.DiGraph]) -> nx.DiGraph:
    if not graphs:
        return nx.DiGraph()

    directed_graphs = [_as_digraph(graph) for graph in graphs if graph is not None]
    if not directed_graphs:
        return nx.DiGraph()

    common_nodes = set(directed_graphs[0].nodes())
    common_edges = set(directed_graphs[0].edges())
    for graph in directed_graphs[1:]:
        common_nodes &= set(graph.nodes())
        common_edges &= set(graph.edges())

    result = nx.DiGraph()
    result.add_nodes_from(common_nodes)
    for u, v in common_edges:
        if u in common_nodes and v in common_nodes:
            result.add_edge(u, v)
    return result


def _reachable_shortest_path_ratio(graph: nx.DiGraph) -> str:
    node_count = graph.number_of_nodes()
    if node_count <= 1:
        return "0 (0%)"

    total_pairs = node_count * (node_count - 1)
    reachable_pairs = 0
    for source in graph.nodes():
        lengths = nx.single_source_shortest_path_length(graph, source)
        reachable_pairs += max(0, len(lengths) - 1)

    ratio = (reachable_pairs / total_pairs) * 100.0
    return f"{reachable_pairs} ({ratio:.0f}%)"


def _compute_network_row(network_name: str, graph: nx.DiGraph) -> SummaryRow:
    node_count = graph.number_of_nodes()
    edge_count = graph.number_of_edges()
    shortest_paths_text = _reachable_shortest_path_ratio(graph)

    if node_count == 0:
        return (network_name, 0, 0, 0, 0, 0.0, shortest_paths_text, 0.0, 0.0)

    largest_component = _largest_undirected_component(graph)
    if largest_component.number_of_nodes() <= 1:
        radius = 0
        diameter = 0
        characteristic_path = 0.0
    else:
        radius = int(nx.radius(largest_component))
        diameter = int(nx.diameter(largest_component))
        characteristic_path = float(nx.average_shortest_path_length(largest_component))

    if node_count > 0:
        avg_neighbors = float(sum(dict(graph.to_undirected().degree()).values()) / node_count)
    else:
        avg_neighbors = 0.0

    clustering_coefficient = float(nx.average_clustering(graph.to_undirected())) if node_count > 0 else 0.0

    return (
        network_name,
        node_count,
        edge_count,
        radius,
        diameter,
        characteristic_path,
        shortest_paths_text,
        avg_neighbors,
        clustering_coefficient,
    )


def build_network_summary_rows(
    initial_graph: nx.Graph | nx.DiGraph,
    pair_graphs: Sequence[nx.Graph | nx.DiGraph],
) -> List[SummaryRow]:
    """Build summary rows for Initial, Union, and Intersection networks."""
    initial_directed = _as_digraph(initial_graph) if initial_graph is not None else nx.DiGraph()
    union_graph = build_union_graph(pair_graphs)
    intersection_graph = build_intersection_graph(pair_graphs)

    return [
        _compute_network_row("Initial", initial_directed),
        _compute_network_row("Union", union_graph),
        _compute_network_row("Intersection", intersection_graph),
    ]
