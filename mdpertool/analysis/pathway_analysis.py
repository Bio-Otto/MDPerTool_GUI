"""Pathway and critical-residue analysis utilities."""

from __future__ import annotations

from typing import Dict, Iterable, List, Sequence, Tuple

import networkx as nx

PathwayRow = Tuple[str, int, object, str]
CriticalRow = Tuple[str, int, float, float]


def summarize_target_pathways(
    source_residue: str,
    targets: Sequence[str],
    graphs: Sequence[nx.Graph],
    progress_callback=None,
) -> Tuple[List[PathwayRow], Dict[str, int], List[str], List[str]]:
    """Build per-target pathway rows and collect intermediate residue usage.

    Returns
    -------
    pathway_rows:
        (target, node_count, shortest_path_length|"N/A", route_type)
    residue_path_hits:
        counts for intermediate residues appearing in shortest paths
    shortest_path_strings:
        human-readable shortest-path strings for list widget
    all_paths_messages:
        detail messages used in DONE summary dialog
    """
    pathway_rows: List[PathwayRow] = []
    residue_path_hits: Dict[str, int] = {}
    shortest_path_strings: List[str] = []
    all_paths_messages: List[str] = []

    total_targets = max(1, len(targets))

    for index, (target_residue, graph) in enumerate(zip(targets, graphs), start=1):
        if source_residue in graph and target_residue in graph and nx.has_path(graph, source_residue, target_residue):
            shortest_path = nx.shortest_path(graph, source_residue, target_residue)
            path_length = max(0, len(shortest_path) - 1)
            route_type = "Direct" if path_length <= 2 else "Indirect"
            pathway_rows.append((target_residue, len(graph.nodes()), path_length, route_type))

            shortest_path_strings.append(" --> ".join(shortest_path))
            all_paths_messages.append(
                "Source: %s  Target: %s\nTotal node number of source-target pair network is : %s"
                % (source_residue, target_residue, len(graph.nodes()))
            )

            for residue_name in shortest_path:
                if residue_name != source_residue and residue_name != target_residue:
                    residue_path_hits[residue_name] = residue_path_hits.get(residue_name, 0) + 1
        else:
            node_count = len(graph.nodes()) if hasattr(graph, "nodes") else 0
            pathway_rows.append((target_residue, node_count, "N/A", "No Path"))

        if callable(progress_callback):
            progress_callback(index, total_targets)

    pathway_rows.sort(key=lambda row: (row[3] == "No Path", str(row[2]), row[0]))
    return pathway_rows, residue_path_hits, shortest_path_strings, all_paths_messages


def build_critical_residue_rows(
    residue_path_hits: Dict[str, int],
    global_betweenness: Dict[str, float],
    top_n: int = 20,
) -> List[CriticalRow]:
    """Rank residues by combined path-hit and betweenness score."""
    if not residue_path_hits:
        return []

    max_hits = max(residue_path_hits.values()) if residue_path_hits else 1
    max_betweenness = max(global_betweenness.values()) if global_betweenness else 1.0

    critical_rows: List[CriticalRow] = []
    for residue_name, hit_count in residue_path_hits.items():
        bw_value = float(global_betweenness.get(residue_name, 0.0))
        hit_score = float(hit_count) / float(max_hits)
        bw_score = float(bw_value) / float(max_betweenness) if max_betweenness > 0 else 0.0
        composite_score = (0.6 * hit_score) + (0.4 * bw_score)
        critical_rows.append((residue_name, hit_count, bw_value, composite_score))

    critical_rows.sort(key=lambda row: (-row[3], -row[1], row[0]))
    return critical_rows[:top_n]


def count_reachability(pathway_rows: Iterable[PathwayRow]) -> Tuple[int, int]:
    """Return (reachable, unreachable) target counts."""
    rows = list(pathway_rows)
    reachable_count = sum(1 for _, _, _, route_type in rows if route_type != "No Path")
    unreachable_count = len(rows) - reachable_count
    return reachable_count, unreachable_count


def extract_target_graph_pairs(
    clean_log_list: Sequence[str],
    all_graph_list: Sequence[nx.Graph],
    amino_acid_residue_codes: Sequence[str],
) -> Tuple[List[str], List[nx.Graph]]:
    """Extract target residue labels and matching graphs from log lines."""
    targets: List[str] = []
    graphs: List[nx.Graph] = []

    for idx, log_line in enumerate(clean_log_list):
        if "Target:" not in log_line:
            continue

        target_residue = log_line.split("Target:")[1].split()[0]
        if target_residue[:3] in amino_acid_residue_codes:
            targets.append(target_residue)
            graphs.append(all_graph_list[idx])

    return targets, graphs


def build_done_message(all_paths_messages: Sequence[str]) -> str:
    """Build DONE message body from per-target lines."""
    return "".join(f"\n{line}" for line in all_paths_messages if line)
