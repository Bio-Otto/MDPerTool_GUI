import os

import numpy as np
import pandas as pd


RESPONSE_THRESHOLD = 0.01


def _build_cumulative_response_curve(response_time_array, num_frames):
    if num_frames <= 0:
        return np.asarray([]), np.asarray([])

    frames = np.arange(num_frames, dtype=float)
    responded_times = response_time_array[response_time_array < num_frames]

    if responded_times.size == 0:
        return frames, np.zeros(num_frames, dtype=float)

    counts_per_frame = np.bincount(responded_times.astype(int), minlength=num_frames)
    cumulative = np.cumsum(counts_per_frame).astype(float)
    return frames, cumulative


def _fit_logistic_cumulative(frames, cumulative, total_residues):
    if total_residues <= 0 or frames.size == 0:
        return {
            "k_d": np.nan,
            "t_half_fit": np.nan,
            "t_half_empirical": np.nan,
            "rmse": np.nan,
            "fitted_curve": np.asarray([]),
        }

    normalized = cumulative / float(total_residues)
    eps = 1e-6
    mask = (normalized > eps) & (normalized < 1.0 - eps)

    t_half_empirical = np.nan
    half_level = total_residues / 2.0
    half_idx = np.where(cumulative >= half_level)[0]
    if half_idx.size > 0:
        t_half_empirical = float(frames[half_idx[0]])

    if np.count_nonzero(mask) < 2:
        return {
            "k_d": np.nan,
            "t_half_fit": np.nan,
            "t_half_empirical": t_half_empirical,
            "rmse": np.nan,
            "fitted_curve": np.full_like(frames, np.nan, dtype=float),
        }

    x = frames[mask]
    y = normalized[mask]
    logit = np.log(y / (1.0 - y))

    slope, intercept = np.polyfit(x, logit, 1)
    if np.isclose(slope, 0.0):
        return {
            "k_d": 0.0,
            "t_half_fit": np.nan,
            "t_half_empirical": t_half_empirical,
            "rmse": np.nan,
            "fitted_curve": np.full_like(frames, np.nan, dtype=float),
        }

    t_half_fit = -intercept / slope
    fitted_norm = 1.0 / (1.0 + np.exp(-(slope * frames + intercept)))
    fitted_curve = fitted_norm * float(total_residues)
    rmse = float(np.sqrt(np.mean((fitted_curve - cumulative) ** 2)))

    return {
        "model": "logistic",
        "k_d": float(slope),
        "t_half_fit": float(t_half_fit),
        "t_half_empirical": t_half_empirical,
        "rmse": rmse,
        "fitted_curve": fitted_curve,
    }


def _fit_gompertz_cumulative(frames, cumulative, total_residues):
    if total_residues <= 0 or frames.size == 0:
        return {
            "model": "gompertz",
            "k_d": np.nan,
            "t_half_fit": np.nan,
            "t_half_empirical": np.nan,
            "rmse": np.nan,
            "fitted_curve": np.asarray([]),
        }

    normalized = cumulative / float(total_residues)
    eps = 1e-6
    mask = (normalized > eps) & (normalized < 1.0 - eps)

    t_half_empirical = np.nan
    half_level = total_residues / 2.0
    half_idx = np.where(cumulative >= half_level)[0]
    if half_idx.size > 0:
        t_half_empirical = float(frames[half_idx[0]])

    if np.count_nonzero(mask) < 2:
        return {
            "model": "gompertz",
            "k_d": np.nan,
            "t_half_fit": np.nan,
            "t_half_empirical": t_half_empirical,
            "rmse": np.nan,
            "fitted_curve": np.full_like(frames, np.nan, dtype=float),
        }

    x = frames[mask]
    y = normalized[mask]

    transformed = np.log(-np.log(y))
    slope, intercept = np.polyfit(x, transformed, 1)

    if np.isclose(slope, 0.0):
        return {
            "model": "gompertz",
            "k_d": 0.0,
            "t_half_fit": np.nan,
            "t_half_empirical": t_half_empirical,
            "rmse": np.nan,
            "fitted_curve": np.full_like(frames, np.nan, dtype=float),
        }

    target = np.log(-np.log(0.5))
    t_half_fit = (target - intercept) / slope

    exponent = np.clip(slope * frames + intercept, -700, 700)
    fitted_norm = np.exp(-np.exp(exponent))
    fitted_curve = fitted_norm * float(total_residues)
    rmse = float(np.sqrt(np.mean((fitted_curve - cumulative) ** 2)))

    return {
        "model": "gompertz",
        "k_d": float(-slope),
        "t_half_fit": float(t_half_fit),
        "t_half_empirical": t_half_empirical,
        "rmse": rmse,
        "fitted_curve": fitted_curve,
    }


def _compute_aic(observed_values, fitted_values, parameter_count=2):
    if observed_values.size == 0 or fitted_values.size == 0:
        return np.nan

    residuals = observed_values - fitted_values
    rss = float(np.sum(residuals ** 2))
    n = int(observed_values.size)

    if n <= 0 or rss <= 0.0:
        return np.nan

    return float(n * np.log(rss / n) + 2.0 * float(parameter_count))


def _select_best_cumulative_fit(frames, cumulative, total_residues):
    logistic_fit = _fit_logistic_cumulative(frames, cumulative, total_residues)
    gompertz_fit = _fit_gompertz_cumulative(frames, cumulative, total_residues)

    aic_logistic = _compute_aic(cumulative, logistic_fit["fitted_curve"], parameter_count=2)
    aic_gompertz = _compute_aic(cumulative, gompertz_fit["fitted_curve"], parameter_count=2)

    candidate_scores = {
        "logistic": aic_logistic,
        "gompertz": aic_gompertz,
    }
    finite_scores = {key: value for key, value in candidate_scores.items() if np.isfinite(value)}

    if finite_scores:
        selected_model = min(finite_scores, key=finite_scores.get)
    else:
        selected_model = "none"

    if selected_model == "logistic":
        selected_fit = logistic_fit
    elif selected_model == "gompertz":
        selected_fit = gompertz_fit
    else:
        selected_fit = {
            "model": "none",
            "k_d": np.nan,
            "t_half_fit": np.nan,
            "t_half_empirical": logistic_fit.get("t_half_empirical", np.nan),
            "rmse": np.nan,
            "fitted_curve": np.full_like(frames, np.nan, dtype=float),
        }

    return {
        "selected_model": selected_model,
        "selected_fit": selected_fit,
        "logistic_fit": logistic_fit,
        "gompertz_fit": gompertz_fit,
        "aic_logistic": aic_logistic,
        "aic_gompertz": aic_gompertz,
    }


def _write_response_metrics(output_name, response_time_array, num_frames):
    frames, cumulative = _build_cumulative_response_curve(response_time_array, num_frames)
    fit_bundle = _select_best_cumulative_fit(frames, cumulative, total_residues=response_time_array.size)
    fit = fit_bundle["selected_fit"]

    responded_count = int(np.count_nonzero(response_time_array < num_frames))
    non_responded_count = int(response_time_array.size - responded_count)
    responded_fraction = float(responded_count / response_time_array.size) if response_time_array.size else 0.0

    na_reason_code = "none"
    fit_status = "ok"
    if responded_count == 0:
        na_reason_code = "all_nonresponsive"
        fit_status = "unavailable"
    elif fit_bundle["selected_model"] == "none":
        na_reason_code = "insufficient_fit_points"
        fit_status = "unavailable"
    elif np.isnan(fit["k_d"]) or np.isnan(fit["t_half_fit"]) or np.isnan(fit["rmse"]):
        na_reason_code = "insufficient_fit_points"
        fit_status = "unavailable"

    summary = pd.DataFrame([
        {
            "total_residues": int(response_time_array.size),
            "responded_residues": responded_count,
            "non_responded_residues": non_responded_count,
            "max_frame": int(max(num_frames - 1, 0)),
            "t_half_empirical_frame": fit["t_half_empirical"],
            "t_half_fit_frame": fit["t_half_fit"],
            "k_d": fit["k_d"],
            "fit_rmse": fit["rmse"],
            "selected_model": fit_bundle["selected_model"],
            "aic_logistic": fit_bundle["aic_logistic"],
            "aic_gompertz": fit_bundle["aic_gompertz"],
            "responded_fraction": responded_fraction,
            "fit_status": fit_status,
            "na_reason_code": na_reason_code,
        }
    ])

    stem, ext = os.path.splitext(output_name)
    summary_path = f"{stem}_metrics.csv"
    fit_curve_path = f"{stem}_fit_curve.csv"
    summary.to_csv(summary_path, index=False)

    fit_curve_df = pd.DataFrame(
        {
            "frame": frames,
            "cumulative_responded_observed": cumulative,
            "cumulative_responded_logistic": fit_bundle["logistic_fit"]["fitted_curve"],
            "cumulative_responded_gompertz": fit_bundle["gompertz_fit"]["fitted_curve"],
            "cumulative_responded_fitted": fit["fitted_curve"],
            "selected_model": fit_bundle["selected_model"],
        }
    )
    fit_curve_df.to_csv(fit_curve_path, index=False)

    return summary_path, fit_curve_path


def get_residue_response_times(reference_name, perturbed_name, output_name='responseTimes.csv'):
    """Calculate the first responsive frame for each residue and write it to disk."""

    reference_energies = pd.read_csv(reference_name)
    perturbed_energies = pd.read_csv(perturbed_name)

    energy_diff = reference_energies.sub(perturbed_energies)
    energy_diff = energy_diff.mask(energy_diff.abs() < RESPONSE_THRESHOLD, 0.0)

    num_frames = len(energy_diff.index)
    residue_response_times = []

    for column_name in energy_diff.columns:
        residue_diff = energy_diff[column_name].to_numpy()
        non_zero_indices = np.flatnonzero(residue_diff)

        if non_zero_indices.size:
            residue_response_times.append(int(non_zero_indices[0]))
        else:
            residue_response_times.append(num_frames)

    response_time_array = np.asarray(residue_response_times)
    np.savetxt(output_name, response_time_array, delimiter=',')
    _write_response_metrics(output_name, response_time_array, num_frames)
    return response_time_array


def getResidueResponseTimes(referenceName, perturbedName, outputName='responseTimes.csv'):
    """Backward-compatible wrapper for legacy callers."""

    return get_residue_response_times(referenceName, perturbedName, outputName)


def summarize_response_time_group(response_time_files, output_name=None):
    """Build a group-level summary across multiple responseTimes CSV files."""

    if not response_time_files:
        raise ValueError("response_time_files cannot be empty")

    filtered_files = []
    for response_file in response_time_files:
        basename = os.path.basename(response_file)
        if basename.endswith("_metrics.csv"):
            continue
        if basename.endswith("_fit_curve.csv"):
            continue
        if basename == "group_summary.csv":
            continue
        filtered_files.append(response_file)

    if not filtered_files:
        raise ValueError("No raw responseTimes CSV files found after filtering sidecars")

    metric_rows = []
    for response_file in filtered_files:
        response_array = np.atleast_1d(np.loadtxt(response_file, delimiter=',')).astype(int)
        num_frames = (int(np.max(response_array)) + 1) if response_array.size else 0

        summary_path, _ = _write_response_metrics(response_file, response_array, num_frames)
        summary_df = pd.read_csv(summary_path)
        if summary_df.empty:
            continue

        row = summary_df.iloc[0].to_dict()
        row["response_time_file"] = os.path.basename(response_file)
        metric_rows.append(row)

    if not metric_rows:
        raise RuntimeError("No valid metrics rows could be generated for group summary")

    metrics_df = pd.DataFrame(metric_rows)

    value_columns = [
        "total_residues",
        "responded_residues",
        "non_responded_residues",
        "max_frame",
        "t_half_empirical_frame",
        "t_half_fit_frame",
        "k_d",
        "fit_rmse",
        "aic_logistic",
        "aic_gompertz",
        "responded_fraction",
    ]

    for column in value_columns:
        if column in metrics_df.columns:
            metrics_df[column] = pd.to_numeric(metrics_df[column], errors='coerce')

    summary_record = {
        "n_files": int(len(metrics_df.index)),
        "fit_ok_count": int((metrics_df.get("fit_status", pd.Series(dtype=str)) == "ok").sum()),
    }

    for column in value_columns:
        if column in metrics_df.columns:
            summary_record[f"{column}_mean"] = float(metrics_df[column].mean(skipna=True))
            summary_record[f"{column}_std"] = float(metrics_df[column].std(skipna=True, ddof=0))

    if "na_reason_code" in metrics_df.columns:
        for reason_code, count in metrics_df["na_reason_code"].value_counts(dropna=False).items():
            key = f"na_reason_{str(reason_code)}_count"
            summary_record[key] = int(count)

    if "selected_model" in metrics_df.columns:
        for model_name, count in metrics_df["selected_model"].value_counts(dropna=False).items():
            key = f"selected_model_{str(model_name)}_count"
            summary_record[key] = int(count)

    if output_name is None:
        output_dir = os.path.dirname(os.path.abspath(filtered_files[0]))
        output_name = os.path.join(output_dir, "group_summary.csv")

    pd.DataFrame([summary_record]).to_csv(output_name, index=False)
    return output_name


__all__ = [
    "RESPONSE_THRESHOLD",
    "get_residue_response_times",
    "getResidueResponseTimes",
    "summarize_response_time_group",
]


# getResidueResponseTimes('reference_energy_file.csv', 'modified_energy_file.csv')
