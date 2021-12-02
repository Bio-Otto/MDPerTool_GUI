from PySide2.QtWidgets import *
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from ast import literal_eval

class WidgetPlot(QWidget):
    def __init__(self, *args, **kwargs):
        QWidget.__init__(self, *args, **kwargs)
        self.canvas = PlotCanvas(self, width=10, height=8)
        self.toolbar = NavigationToolbar(self.canvas, self)


class PlotCanvas(FigureCanvas):
    def __init__(self, parent=None, width=10, height=8, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        FigureCanvas.__init__(self, fig)
        self.setParent(parent)
        FigureCanvas.setSizePolicy(self, QSizePolicy.Expanding, QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)
        # self.plot()

    def plot(self, data, source_residue=None):
        ax = self.figure.add_subplot(111)
        ax.plot(data, 'r-', linewidth=0.5)
        ax.set_title('Response Time Graph')

        ax.set_xlabel('Time Step (Frame)')
        ax.set_ylabel('Responded Residue Count')

        max_handler = 0
        for i, line in enumerate(ax.lines):
            x_data = ax.lines[i].get_xdata()  # get the first line, there might be more
            try:
                max_val = max(x_data)
                if max_val > max_handler:
                    max_handler = max_val

            except Exception as ERr:
                print(ERr)

        if source_residue is not None:
            ax.text(max_handler - 350, 10, 'Perturbed Residue(s): %s' % source_residue, style='italic',
                     bbox={'facecolor': 'red', 'alpha': 0.5, 'pad': 5})
            ax.legend(title='Speed Factor')

        ax.margins(x=0.01, y=0.01, tight=True)
        self.figure.tight_layout()


        # fig.set_size_inches((25, 15), forward=False)
        # fig.savefig('%s.png' % top_file_name, dpi=300)  # Change is over here

        self.draw()


def getResponseTimeGraph(responseTimeFile):
    # Draw Graph residue response vs Time
    Responses_file = pd.read_csv(responseTimeFile, header=None)
    Response_Time_Column = Responses_file.values
    Response_Time_Column_Max = Response_Time_Column.max()
    row, col = Responses_file.shape
    Response_Count = []
    Increase = 0

    for i in range(int(Response_Time_Column_Max)):
        Increase = Increase + list(Response_Time_Column).count(i)
        Response_Count.append(Increase)
    # plt.ylim(0, row + 5)
    # plt.plot(Response_Count, linewidth=1)

    return row, col, Response_Count

    """
    if speed_factors is not None:
        for factor in sorted(speed_factors):
            try:
                row_num, col_num = getResponseTimeGraph('responseTimes_%s.csv' % factor, factor=factor)
            except:
                print("There is no % named file in the directory")


    """