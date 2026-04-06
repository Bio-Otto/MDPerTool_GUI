"""Service helpers for preparing response-dynamics analysis outputs."""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any, Dict, List, Tuple, Mapping, Optional

import numpy as np

from analysis.residue_response_analyzer import ResidueResponseAnalyzer


def build_response_dynamics_payload(
    possible_path: str,
    metrics: Mapping[Any, Any],
    frame_time_delta: float = 1.0,
    residue_names: Optional[List[str]] = None,
) -> Dict[str, Any]:
    """Build UI-friendly payloads for response-dynamics tables."""
    payload: Dict[str, Any] = {
        "residue_rows": [],
        "domain_rows": [],
        "metrics_rows": [],
        "qc_rows": [],
        "provenance_rows": [],
    }

    response_path = Path(possible_path)
    response_stem = response_path.with_suffix("")
    response_metrics_file = str(response_stem) + "_metrics.csv"
    response_fit_file = str(response_stem) + "_fit_curve.csv"

    if os.path.isfile(possible_path) and os.path.isfile(response_metrics_file) and os.path.isfile(response_fit_file):
        response_analyzer = ResidueResponseAnalyzer(
            possible_path,
            response_metrics_file,
            response_fit_file,
            frame_time_delta=frame_time_delta,
        )

        if residue_names is None:
            residue_names = [f"RES_{i + 1}" for i in range(response_analyzer.num_residues)]

        residue_summary_df = response_analyzer.get_per_residue_summary(residue_names)
        payload["residue_rows"] = [
            (
                str(row_data["residue_id"]),
                str(row_data["residue_name"]),
                str(row_data["response_time_frame"]),
            )
            for _, row_data in residue_summary_df.iterrows()
        ]

        groups = response_analyzer.get_residue_groups_by_threshold()
        payload["domain_rows"] = [
            (
                "Fast Responders",
                len(groups["fast"]),
                np.mean(response_analyzer.response_times_ps[groups["fast"]]) if groups["fast"] else 0.0,
                np.min(response_analyzer.response_times_ps[groups["fast"]]) if groups["fast"] else 0.0,
                np.max(response_analyzer.response_times_ps[groups["fast"]]) if groups["fast"] else 0.0,
            ),
            (
                "Medium Responders",
                len(groups["medium"]),
                np.mean(response_analyzer.response_times_ps[groups["medium"]]) if groups["medium"] else 0.0,
                np.min(response_analyzer.response_times_ps[groups["medium"]]) if groups["medium"] else 0.0,
                np.max(response_analyzer.response_times_ps[groups["medium"]]) if groups["medium"] else 0.0,
            ),
            (
                "Slow Responders",
                len(groups["slow"]),
                np.mean(response_analyzer.response_times_ps[groups["slow"]]) if groups["slow"] else 0.0,
                np.min(response_analyzer.response_times_ps[groups["slow"]]) if groups["slow"] else 0.0,
                np.max(response_analyzer.response_times_ps[groups["slow"]]) if groups["slow"] else 0.0,
            ),
        ]

    metrics_order = [
        ("total_residues", "Total Residues"),
        ("responded_residues", "Responded Residues"),
        ("non_responded_residues", "Non-Responded Residues"),
        ("max_frame", "Max Frame"),
        ("t_half_empirical_frame", "t_half Empirical"),
        ("t_half_fit_frame", "t_half Fit"),
        ("k_d", "k_d"),
        ("fit_rmse", "Fit RMSE"),
        ("selected_model", "Selected Model"),
        ("aic_logistic", "AIC Logistic"),
        ("aic_gompertz", "AIC Gompertz"),
    ]

    metrics_rows: List[Tuple[str, str]] = []
    for metric_key, metric_label in metrics_order:
        value = metrics.get(metric_key)
        if isinstance(value, float):
            value_text = f"{value:.4f}"
        elif value is None:
            value_text = "N/A"
        else:
            value_text = str(value)
        metrics_rows.append((metric_label, value_text))
    payload["metrics_rows"] = metrics_rows

    responded_fraction = metrics.get("responded_fraction")
    fit_status = str(metrics.get("fit_status") or "unavailable")
    na_reason_code = str(metrics.get("na_reason_code") or "none")
    fit_rmse = metrics.get("fit_rmse")

    if isinstance(fit_rmse, float):
        if fit_rmse < 10.0:
            rmse_level = "low"
        elif fit_rmse < 30.0:
            rmse_level = "medium"
        else:
            rmse_level = "high"
    else:
        rmse_level = "N/A"

    if isinstance(responded_fraction, float):
        responded_fraction_text = f"{responded_fraction * 100.0:.1f}%"
    else:
        responded_fraction_text = "N/A"

    payload["qc_rows"] = [
        ("fit_ok", "yes" if fit_status == "ok" else "no"),
        ("rmse_level", rmse_level),
        ("responded_fraction", responded_fraction_text),
        ("na_reason_code", na_reason_code),
    ]

    manifest_path = f"{os.path.splitext(possible_path)[0]}_analysis_manifest.json"
    provenance_data = {
        "manifest_exists": "no",
        "created_at_utc": "N/A",
        "selected_model": str(metrics.get("selected_model") or "N/A"),
        "fit_status": fit_status,
    }

    if os.path.exists(manifest_path):
        try:
            with open(manifest_path, "r", encoding="utf-8") as manifest_file:
                manifest_json = json.load(manifest_file)
            provenance_data["manifest_exists"] = "yes"
            provenance_data["created_at_utc"] = str(manifest_json.get("created_at_utc") or "N/A")
            provenance_data["selected_model"] = str(
                manifest_json.get("selected_model") or provenance_data["selected_model"]
            )
            provenance_data["fit_status"] = str(
                manifest_json.get("fit_status") or provenance_data["fit_status"]
            )
        except Exception:
            provenance_data["manifest_exists"] = "error"

    payload["provenance_rows"] = [
        ("manifest_exists", provenance_data["manifest_exists"]),
        ("created_at_utc", provenance_data["created_at_utc"]),
        ("selected_model", provenance_data["selected_model"]),
        ("manifest_path", manifest_path),
    ]

    return payload
