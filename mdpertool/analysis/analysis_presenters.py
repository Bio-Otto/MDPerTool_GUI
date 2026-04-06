"""UI presenter helpers for analysis tables and warnings."""

from __future__ import annotations

from typing import Any, Mapping, Sequence, Tuple

from PySide2 import QtCore
from PySide2.QtGui import QColor, QBrush
from PySide2.QtWidgets import QTableWidget, QTableWidgetItem


def present_warning_payload(
    show_warning_fn: Any,
    owner: Any,
    payload: Mapping[str, Any],
    stylesheet: Any,
    default_title: str = "Warning",
    default_message: str = "Unexpected warning.",
) -> None:
    """Display a warning payload via provided warning presenter function."""
    title = str(payload.get("title", default_title))
    message = str(payload.get("message", default_message))
    show_warning_fn(owner, title, message, stylesheet)


def populate_pathway_summary_table(table: QTableWidget, rows: Sequence[Tuple[str, int, object, str]]) -> None:
    table.setRowCount(len(rows))
    for row_idx, (target_name, node_count, path_len, route_type) in enumerate(rows):
        table.setItem(row_idx, 0, QTableWidgetItem(str(target_name)))
        table.setItem(row_idx, 1, QTableWidgetItem(str(node_count)))
        table.setItem(row_idx, 2, QTableWidgetItem(str(path_len)))
        table.setItem(row_idx, 3, QTableWidgetItem(route_type))


def populate_critical_residue_table(table: QTableWidget, rows: Sequence[Tuple[str, int, float, float]]) -> None:
    table.setRowCount(len(rows))
    for row_idx, (residue_name, hit_count, bw_value, composite_score) in enumerate(rows):
        table.setItem(row_idx, 0, QTableWidgetItem(str(residue_name)))
        table.setItem(row_idx, 1, QTableWidgetItem(str(hit_count)))
        table.setItem(row_idx, 2, QTableWidgetItem(f"{bw_value:.6f}"))
        table.setItem(row_idx, 3, QTableWidgetItem(f"{composite_score:.4f}"))


def populate_reachability_qc(table: QTableWidget, reachable_count: int, unreachable_count: int, start_row: int = 4) -> None:
    table.setItem(start_row, 0, QTableWidgetItem("reachable_targets"))
    table.setItem(start_row, 1, QTableWidgetItem(str(reachable_count)))
    table.setItem(start_row + 1, 0, QTableWidgetItem("unreachable_targets"))
    table.setItem(start_row + 1, 1, QTableWidgetItem(str(unreachable_count)))


def populate_residue_response_table(table: QTableWidget, rows: Sequence[Tuple[str, str, str]]) -> None:
    table.setRowCount(len(rows))
    for row_idx, (residue_id, residue_name, response_frame) in enumerate(rows):
        table.setItem(row_idx, 0, QTableWidgetItem(str(residue_id)))
        table.setItem(row_idx, 1, QTableWidgetItem(str(residue_name)))
        table.setItem(row_idx, 2, QTableWidgetItem(str(response_frame)))


def populate_domain_summary_table(table: QTableWidget, rows: Sequence[Tuple[str, int, float, float, float]]) -> None:
    table.setRowCount(len(rows))
    for row_idx, (domain_name, count, mean_time, min_time, max_time) in enumerate(rows):
        table.setItem(row_idx, 0, QTableWidgetItem(str(domain_name)))
        table.setItem(row_idx, 1, QTableWidgetItem(str(count)))
        table.setItem(row_idx, 2, QTableWidgetItem(f"{mean_time:.2f}"))
        table.setItem(row_idx, 3, QTableWidgetItem(f"{min_time:.2f}"))
        table.setItem(row_idx, 4, QTableWidgetItem(f"{max_time:.2f}"))


def populate_network_summary_table(
    table: QTableWidget,
    rows: Sequence[Tuple[str, int, int, int, int, float, str, float, float]],
) -> None:
    table.setRowCount(len(rows))
    for row_idx, row in enumerate(rows):
        (
            network_name,
            node_count,
            edge_count,
            radius,
            diameter,
            characteristic_path_length,
            shortest_paths_text,
            average_neighbors,
            clustering_coefficient,
        ) = row

        table.setItem(row_idx, 0, QTableWidgetItem(str(network_name)))
        table.setItem(row_idx, 1, QTableWidgetItem(str(node_count)))
        table.setItem(row_idx, 2, QTableWidgetItem(str(edge_count)))
        table.setItem(row_idx, 3, QTableWidgetItem(str(radius)))
        table.setItem(row_idx, 4, QTableWidgetItem(str(diameter)))
        table.setItem(row_idx, 5, QTableWidgetItem(f"{characteristic_path_length:.3f}"))
        table.setItem(row_idx, 6, QTableWidgetItem(str(shortest_paths_text)))
        table.setItem(row_idx, 7, QTableWidgetItem(f"{average_neighbors:.3f}"))
        table.setItem(row_idx, 8, QTableWidgetItem(f"{clustering_coefficient:.3f}"))


def populate_motif_summary_table(
    table: QTableWidget,
    rows: Sequence[Tuple[str, str, int, int, str, str]],
) -> None:
    table.setRowCount(len(rows))
    for row_idx, (size_name, motif_id, edge_count, occurrence, frequency_pct, scope_name) in enumerate(rows):
        table.setItem(row_idx, 0, QTableWidgetItem(str(size_name)))
        table.setItem(row_idx, 1, QTableWidgetItem(str(motif_id)))
        table.setItem(row_idx, 2, QTableWidgetItem(str(edge_count)))
        table.setItem(row_idx, 3, QTableWidgetItem(str(occurrence)))
        table.setItem(row_idx, 4, QTableWidgetItem(str(frequency_pct)))
        table.setItem(row_idx, 5, QTableWidgetItem(str(scope_name)))


def populate_metrics_table(table: QTableWidget, rows: Sequence[Tuple[str, str]]) -> None:
    table.setRowCount(len(rows))
    for row_idx, (metric_label, value_text) in enumerate(rows):
        metric_item = QTableWidgetItem(str(metric_label))
        value_item = QTableWidgetItem(str(value_text))
        metric_item.setFlags(QtCore.Qt.ItemIsEnabled)
        value_item.setFlags(QtCore.Qt.ItemIsEnabled)
        table.setItem(row_idx, 0, metric_item)
        table.setItem(row_idx, 1, value_item)


def populate_qc_table(table: QTableWidget, rows: Sequence[Tuple[str, str]]) -> None:
    table.setRowCount(max(table.rowCount(), len(rows)))
    reason_text_map = {
        "none": "No NA condition detected.",
        "missing_sidecar": "Metrics sidecar was missing or could not be loaded.",
        "insufficient_fit_points": "Insufficient sigmoid-transition points for fit.",
        "all_nonresponsive": "No residue crossed response threshold.",
    }

    for row_idx, (qc_key, qc_value) in enumerate(rows):
        qc_key_item = QTableWidgetItem(str(qc_key))
        qc_value_item = QTableWidgetItem(str(qc_value))
        qc_key_item.setFlags(QtCore.Qt.ItemIsEnabled)
        qc_value_item.setFlags(QtCore.Qt.ItemIsEnabled)

        if qc_key == "fit_ok":
            if str(qc_value).lower() == "yes":
                qc_value_item.setBackground(QBrush(QColor("#77dd77")))
                qc_value_item.setToolTip("Fit quality check passed.")
            else:
                qc_value_item.setBackground(QBrush(QColor("#ff686b")))
                qc_value_item.setToolTip("Fit quality check failed or unavailable.")

        if qc_key == "rmse_level":
            rmse_level_value = str(qc_value).lower()
            if rmse_level_value == "low":
                qc_value_item.setBackground(QBrush(QColor("#77dd77")))
                qc_value_item.setToolTip("Low fit error.")
            elif rmse_level_value == "medium":
                qc_value_item.setBackground(QBrush(QColor("#f6bc66")))
                qc_value_item.setToolTip("Moderate fit error.")
            elif rmse_level_value == "high":
                qc_value_item.setBackground(QBrush(QColor("#ff686b")))
                qc_value_item.setToolTip("High fit error; inspect curve carefully.")
            else:
                qc_value_item.setToolTip("RMSE not available.")

        if qc_key == "responded_fraction":
            qc_value_item.setToolTip("Fraction of residues that responded before max frame.")

        if qc_key == "na_reason_code":
            reason_key = str(qc_value)
            qc_value_item.setToolTip(reason_text_map.get(reason_key, "Unknown NA reason code."))

        table.setItem(row_idx, 0, qc_key_item)
        table.setItem(row_idx, 1, qc_value_item)


def populate_provenance_table(table: QTableWidget, rows: Sequence[Tuple[str, str]]) -> None:
    table.setRowCount(len(rows))
    for row_idx, (prov_key, prov_value) in enumerate(rows):
        prov_key_item = QTableWidgetItem(str(prov_key))
        prov_value_item = QTableWidgetItem(str(prov_value))
        prov_key_item.setFlags(QtCore.Qt.ItemIsEnabled)
        prov_value_item.setFlags(QtCore.Qt.ItemIsEnabled)

        if prov_key == "manifest_exists":
            prov_value_lower = str(prov_value).lower()
            if prov_value_lower == "yes":
                prov_value_item.setBackground(QBrush(QColor("#77dd77")))
            elif prov_value_lower == "error":
                prov_value_item.setBackground(QBrush(QColor("#ff686b")))
            else:
                prov_value_item.setBackground(QBrush(QColor("#f6bc66")))

        if prov_key == "manifest_path":
            prov_value_item.setToolTip("Analysis manifest location on disk.")

        table.setItem(row_idx, 0, prov_key_item)
        table.setItem(row_idx, 1, prov_value_item)
