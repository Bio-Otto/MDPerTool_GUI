from PySide2.QtWidgets import QWidget, QSizePolicy
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
import pandas as pd
import os
import warnings
import json
from datetime import datetime, timezone
from Bio import BiopythonWarning
# Suppress specific warnings
warnings.filterwarnings("ignore", category=UserWarning, message="Tight layout not applied.*")
warnings.filterwarnings("ignore", category=UserWarning, message=".*PDBConstructionWarning.*")
warnings.filterwarnings("ignore", category=BiopythonWarning, message=".*PDBConstructionWarning.*")

try:
    from _version import __version__ as mdpertool_version
except Exception:
    try:
        from mdpertool._version import __version__ as mdpertool_version
    except Exception:
        mdpertool_version = "unknown"


class WidgetPlot(QWidget):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.canvas = PlotCanvas(self, width=10, height=5)
        self.toolbar = NavigationToolbar(self.canvas, self)


class PlotCanvas(FigureCanvas):
    def __init__(self, parent=None, width=10, height=5, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        super().__init__(fig)
        self.setParent(parent)
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.updateGeometry()

    def plot(self, data, source_residue=None, plot_name=None, fitted_data=None):
        self.figure.clear()
        ax = self.figure.add_subplot(111)
        ax.plot(data, 'r-', linewidth=1.0, label=plot_name or 'Observed')

        if fitted_data is not None and len(fitted_data):
            ax.plot(fitted_data, color='dodgerblue', linestyle='--', linewidth=1.0, label='Fitted')

        ax.set_title('Response Time Graph')
        ax.set_xlabel('Time Step (Frame)')
        ax.set_ylabel('Responded Residue Count')

        max_val = max(ax.lines[0].get_xdata(), default=0) if ax.lines else 0

        if source_residue is not None:
            ax.text(max(max_val - 470, 0), 18, f'Perturbed Residue(s): {source_residue}', style='italic',
                    bbox={'facecolor': 'red', 'alpha': 0.5, 'pad': 5}, fontsize=8)

        ax.legend(title='Curve')

        ax.margins(x=0.01, y=0.01, tight=True)

        # Adjust the tight_layout to avoid the warning
        try:
            self.figure.tight_layout(rect=(0.1, 0.1, 0.9, 0.9))
        except (RuntimeError, ValueError) as e:
            print(f"Warning: {e}")

        self.draw()


def getResponseTimeGraph(responseTimeFile):
    responses_file = pd.read_csv(responseTimeFile, header=None)
    response_time_column = responses_file.values.flatten()
    response_time_column_max = response_time_column.max()

    response_count = [sum(response_time_column <= i) for i in range(int(response_time_column_max) + 1)]

    row, col = responses_file.shape
    filename_parts = os.path.splitext(responseTimeFile)  # Split filename and extension
    plot_name = "x%s" % str(filename_parts[0].split('_')[-1])  # Extract last segment after last underscore

    fit_curve = []
    metrics = {}
    generation_error = None
    fit_curve_file = f"{filename_parts[0]}_fit_curve.csv"
    metrics_file = f"{filename_parts[0]}_metrics.csv"

    if not (os.path.exists(fit_curve_file) and os.path.exists(metrics_file)):
        try:
            try:
                from no_gui.response_time_creator import _write_response_metrics
            except ModuleNotFoundError:
                from mdpertool.no_gui.response_time_creator import _write_response_metrics

            num_frames = int(response_time_column_max) + 1
            _write_response_metrics(
                responseTimeFile,
                response_time_column.astype(int),
                num_frames,
            )
        except Exception as err:
            generation_error = err
            print(f"Warning: could not auto-generate response metrics sidecar files: {err}")

    if os.path.exists(fit_curve_file):
        fit_curve_df = pd.read_csv(fit_curve_file)
        if 'cumulative_responded_fitted' in fit_curve_df.columns:
            fit_curve_series = fit_curve_df['cumulative_responded_fitted'].dropna()
            fit_curve = fit_curve_series.tolist()

    if os.path.exists(metrics_file):
        metrics_df = pd.read_csv(metrics_file)
        if not metrics_df.empty:
            metric_row = metrics_df.iloc[0].to_dict()
            metrics = {
                key: (None if pd.isna(value) else value)
                for key, value in metric_row.items()
            }

    if not metrics:
        metrics = {
            "na_reason_code": "missing_sidecar",
            "fit_status": "unavailable",
            "responded_fraction": None,
            "fit_rmse": None,
        }

    if generation_error is not None:
        metrics["na_reason_code"] = "missing_sidecar"

    try:
        _write_analysis_manifest(responseTimeFile, metrics)
    except Exception as err:
        print(f"Warning: could not write analysis manifest: {err}")

    return row, col, response_count, plot_name, fit_curve, metrics


def _write_analysis_manifest(response_time_file, metrics):
    stem, _ = os.path.splitext(response_time_file)
    manifest_path = f"{stem}_analysis_manifest.json"

    manifest = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "mdpertool_version": mdpertool_version,
        "response_time_file": response_time_file,
        "metrics_file": f"{stem}_metrics.csv",
        "fit_curve_file": f"{stem}_fit_curve.csv",
        "selected_model": metrics.get("selected_model"),
        "fit_status": metrics.get("fit_status"),
        "na_reason_code": metrics.get("na_reason_code"),
    }

    with open(manifest_path, 'w', encoding='utf-8') as manifest_file:
        json.dump(manifest, manifest_file, indent=2, ensure_ascii=False)




