from PySide2.QtWidgets import QWidget, QSizePolicy
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
import pandas as pd
import os
import warnings
from Bio import BiopythonWarning
# Suppress specific warnings
warnings.filterwarnings("ignore", category=UserWarning, message="Tight layout not applied.*")
warnings.filterwarnings("ignore", category=UserWarning, message=".*PDBConstructionWarning.*")
warnings.filterwarnings("ignore", category=BiopythonWarning, message=".*PDBConstructionWarning.*")


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

    def plot(self, data, source_residue=None, plot_name=None):
        self.figure.clear()
        ax = self.figure.add_subplot(111)
        ax.plot(data, 'r-', linewidth=0.5, label=plot_name)
        ax.set_title('Response Time Graph')
        ax.set_xlabel('Time Step (Frame)')
        ax.set_ylabel('Responded Residue Count')

        max_val = max(ax.lines[0].get_xdata(), default=0) if ax.lines else 0

        if source_residue is not None:
            ax.text(max(max_val - 470, 0), 18, f'Perturbed Residue(s): {source_residue}', style='italic',
                    bbox={'facecolor': 'red', 'alpha': 0.5, 'pad': 5}, fontsize=8)
            ax.legend(title='Speed Factor')

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
    return row, col, response_count, plot_name




