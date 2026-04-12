"""Lightweight motif visualization helpers for network analysis."""

from __future__ import annotations

from collections import OrderedDict
from typing import Any, List, Sequence, Tuple

import networkx as nx

from .network_motif_service import _count_connected_induced_motifs


MotifVisualRow = Tuple[str, str, str, int, int, str, str, str]


def build_motif_visual_rows(
    graph: nx.DiGraph,
    scope_name: str,
    motif_sizes: Sequence[int] = (3, 4),
    max_combinations: int = 400000,
    top_k_per_size: int = 2,
) -> List[MotifVisualRow]:
    """Return the top motif topologies together with their canonical keys."""
    if graph is None:
        return []

    rows: List[MotifVisualRow] = []

    for motif_size in motif_sizes:
        motif_counter, motif_examples, connected_count, status = _count_connected_induced_motifs(
            graph=graph,
            motif_size=motif_size,
            max_combinations=max_combinations,
        )

        if status != "ok" or connected_count == 0:
            continue

        sorted_items = sorted(motif_counter.items(), key=lambda item: item[1], reverse=True)[:top_k_per_size]
        for idx, (motif_key, occurrence) in enumerate(sorted_items, start=1):
            edge_count = motif_key.count("1")
            frequency_pct = (occurrence / connected_count) * 100.0
            example_residues = ", ".join(motif_examples.get(motif_key, ())) or "N/A"
            rows.append(
                (
                    f"{motif_size}-node",
                    f"M{motif_size}_{idx:02d}",
                    motif_key,
                    edge_count,
                    occurrence,
                    f"{frequency_pct:.2f}%",
                    scope_name,
                    example_residues,
                )
            )

    return rows


def _motif_key_to_graph(motif_key: str, motif_size: int) -> nx.DiGraph:
    graph = nx.DiGraph()
    nodes = [str(index + 1) for index in range(motif_size)]
    graph.add_nodes_from(nodes)

    normalized_key = str(motif_key or "")
    expected_length = motif_size * motif_size
    if len(normalized_key) < expected_length:
        normalized_key = normalized_key.ljust(expected_length, "0")

    for row_index in range(motif_size):
        for col_index in range(motif_size):
            if row_index == col_index:
                continue
            key_index = row_index * motif_size + col_index
            if normalized_key[key_index] == "1":
                graph.add_edge(nodes[row_index], nodes[col_index])

    return graph


def render_motif_gallery(widget: Any, rows: Sequence[MotifVisualRow]) -> None:
    """Render a small gallery of top motif topologies into a Matplotlib widget."""
    figure = widget.canvas.figure
    figure.clear()

    if not rows:
        axis = figure.add_subplot(111)
        axis.axis("off")
        axis.text(
            0.5,
            0.5,
            "No motif data available for visualization.",
            ha="center",
            va="center",
            fontsize=11,
            fontweight="bold",
        )
        figure.tight_layout()
        widget.canvas.draw()
        return

    grouped_rows: OrderedDict[str, List[MotifVisualRow]] = OrderedDict()
    for row in rows:
        grouped_rows.setdefault(row[0], []).append(row)

    total_rows = len(grouped_rows)
    max_columns = max(len(items) for items in grouped_rows.values())
    axes = figure.subplots(total_rows, max_columns, squeeze=False)

    for row_index, (size_name, size_rows) in enumerate(grouped_rows.items()):
        for col_index in range(max_columns):
            axis = axes[row_index][col_index]
            axis.set_axis_off()

            if col_index >= len(size_rows):
                continue

            _, motif_id, motif_key, edge_count, occurrence, frequency_text, scope_name, example_residues = size_rows[col_index]
            motif_size = int(size_name.split("-")[0])
            motif_graph = _motif_key_to_graph(motif_key, motif_size)
            positions = nx.circular_layout(motif_graph)

            axis.set_title(
                f"{motif_id} | {occurrence} hits | {frequency_text}",
                fontsize=9,
                fontweight="bold",
            )
            nx.draw_networkx_nodes(
                motif_graph,
                positions,
                ax=axis,
                node_color="#7aa6c2",
                node_size=650,
                edgecolors="#1f1f1f",
                linewidths=0.8,
            )
            nx.draw_networkx_labels(
                motif_graph,
                positions,
                ax=axis,
                font_size=8,
                font_weight="bold",
                font_color="#111111",
            )
            nx.draw_networkx_edges(
                motif_graph,
                positions,
                ax=axis,
                arrows=True,
                arrowstyle="-|>",
                arrowsize=14,
                width=1.2,
                edge_color="#4f4f4f",
                connectionstyle="arc3,rad=0.08",
            )

            axis.text(
                0.5,
                -0.17,
                f"{size_name} motif in {scope_name}\nResidues: {example_residues}",
                transform=axis.transAxes,
                ha="center",
                va="top",
                fontsize=8,
                color="#555555",
            )

    figure.suptitle("Top Motif Topologies", fontsize=12, fontweight="bold")
    figure.tight_layout(rect=(0, 0.02, 1, 0.95))
    widget.canvas.draw()
