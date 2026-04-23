"""Intermolecular propagation coefficient helpers for inter-subunit signal flow."""

from __future__ import annotations

import os
import re
from typing import Dict, List, Sequence, Tuple

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
from threading import Lock


_IPC_PLOT_LOCK = Lock()


IPCRow = Tuple[str, str, int, int, float, str, str]


def _node_chain_label(graph: nx.DiGraph, node: str) -> str:
	"""Return the best available chain identifier for a node label."""
	if node in graph:
		node_data = graph.nodes[node]
		for key in ("chain", "chain_id", "subunit"):
			value = node_data.get(key)
			if value is not None and str(value).strip():
				return str(value).strip()

	label = str(node).strip()
	match = re.match(r"^[A-Za-z]{3}\d+([A-Za-z])$", label)
	if match:
		return match.group(1)
	return ""


def _collect_interface_nodes(graph: nx.DiGraph) -> List[str]:
	"""Return nodes that participate in at least one inter-chain edge."""
	interface_nodes = set()
	for source_node, target_node in graph.edges():
		source_chain = _node_chain_label(graph, source_node)
		target_chain = _node_chain_label(graph, target_node)
		if source_chain and target_chain and source_chain != target_chain:
			interface_nodes.add(str(source_node))
			interface_nodes.add(str(target_node))

	if interface_nodes:
		return sorted(interface_nodes)

	return sorted(str(node) for node in graph.nodes())


def _normalized_ipc_scores(graph: nx.DiGraph, candidate_nodes: Sequence[str]) -> Dict[str, float]:
	raw_scores: Dict[str, float] = {}
	for node in candidate_nodes:
		in_degree = int(graph.in_degree(node))
		out_degree = int(graph.out_degree(node))
		raw_scores[str(node)] = float(in_degree * out_degree)

	total_score = sum(raw_scores.values())
	if total_score <= 0:
		return {node: 0.0 for node in raw_scores}

	return {node: score / total_score for node, score in raw_scores.items()}


def build_intermolecular_propagation_rows(
	graph: nx.DiGraph,
	top_n: int = 20,
) -> List[IPCRow]:
	"""Rank interface residues by their intermolecular propagation coefficient.

	The score is normalized from the node's in-degree / out-degree product and is
	intended for inter-subunit networks where chain transitions matter more than
	global layer bypassing.
	"""
	if graph.number_of_nodes() == 0:
		return []

	candidate_nodes = _collect_interface_nodes(graph)
	ipc_scores = _normalized_ipc_scores(graph, candidate_nodes)
	rows: List[IPCRow] = []

	for node in candidate_nodes:
		in_degree = int(graph.in_degree(node))
		out_degree = int(graph.out_degree(node))
		chain_label = _node_chain_label(graph, node) or "N/A"
		ipc_value = float(ipc_scores.get(str(node), 0.0))

		if in_degree == 0 and out_degree > 0:
			role = "Initiator / Source"
			rationale = "Starts the inter-subunit signal flow."
		elif out_degree == 0 and in_degree > 0:
			role = "Sink / Target"
			rationale = "Receives inter-subunit signal input."
		elif out_degree > in_degree:
			role = "Exporter / Relay"
			rationale = "Forwards signal across the interface."
		elif in_degree > out_degree:
			role = "Receiver / Funnel"
			rationale = "Collects incoming interface signals."
		else:
			role = "Bridge / Relay"
			rationale = "Balances incoming and outgoing interface flow."

		rows.append((str(node), chain_label, in_degree, out_degree, round(ipc_value, 4), role, rationale))

	rows.sort(key=lambda row: (-row[4], -row[2], -row[3], row[0]))
	return rows[:top_n]


def _write_ipc_bar_plot(rows: Sequence[IPCRow], output_image: str) -> None:
	"""Write a compact IPC ranking plot for quick review."""
	if not rows:
		return

	with _IPC_PLOT_LOCK:
		labels = [row[0] for row in rows]
		values = [row[4] for row in rows]

		fig, ax = plt.subplots(figsize=(10, 6))
		ax.bar(labels, values, color="#4c78a8", edgecolor="black")
		ax.tick_params(axis="x", rotation=45, labelsize=9)
		ax.set_ylabel("Intermolecular Propagation Coefficient (IPC)", fontsize=11)
		ax.set_xlabel("Residue", fontsize=11)
		ax.set_title("Top Inter-Subunit Propagation Residues", fontsize=13, fontweight="bold")
		fig.tight_layout()
		fig.savefig(output_image, dpi=300)
		plt.close(fig)


def execute_intermolecular_propagation_workflow(
	graph: nx.DiGraph,
	output_directory: str,
	prefix: str = "",
	top_n: int = 20,
):
	"""Build IPC CSV and plot outputs for the given network."""
	rows = build_intermolecular_propagation_rows(graph=graph, top_n=top_n)
	if not rows:
		return pd.DataFrame(columns=["Residue_ID", "Chain", "In_Degree", "Out_Degree", "IPC", "Assigned_Role", "Biological_Interpretation"])

	os.makedirs(output_directory, exist_ok=True)
	csv_path = os.path.join(output_directory, f"{prefix}intermolecular_propagation_coefficients.csv")
	plot_path = os.path.join(output_directory, f"{prefix}intermolecular_propagation_plot.png")

	df_rows = pd.DataFrame(
		rows,
		columns=[
			"Residue_ID",
			"Chain",
			"In_Degree",
			"Out_Degree",
			"IPC",
			"Assigned_Role",
			"Biological_Interpretation",
		],
	)
	df_rows.to_csv(csv_path, index=False)
	_write_ipc_bar_plot(rows, plot_path)
	return df_rows