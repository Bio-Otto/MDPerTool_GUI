import os
import queue
import threading
import time
import tokenize
from io import StringIO
import numpy as np
import pyqtgraph as pg
from PySide2 import QtCore
from PySide2.QtWidgets import *
from openmm.app import StateDataReporter


# import multiprocessing
# import pystache
# from PySide2 import QtWidgets, QtCore, QtGui
# from pyqtgraph import PlotWidget, plot, dockarea, ProgressDialog
# from PySide2.QtCore import QTimer, QDateTime, Slot
# import subprocess
# import continuous_threading
# from PySide2.QtGui import (QBrush, QColor, QConicalGradient, QCursor, QFont, QFontDatabase, QIcon, QKeySequence,
#                            QLinearGradient, QPalette, QPainter, QPixmap, QRadialGradient)


# ################################################################################################################### #
# #################################################### FUNCTIONS #################################################### #
# ################################################################################################################### #
from typing import List


def queue_reporter_factory(queue):
    """Factory function that returns a dynamically defined OpenMM
    reporter class which reports by sending dicts down a synchronous queue"""

    class QueueStateDataReporter(StateDataReporter):
        """Subclass of StateDataReporter sends its results down a synchronous
        Queue, as opposed to printing them to a file-like object
        """

        def __init__(self, file, *args, **kwargs):
            with open(os.devnull, 'w') as f:
                # send in a fake file
                super(QueueStateDataReporter, self).__init__(f, *args, **kwargs)

            # this is where we'll store the names of the fields that
            # are being reported in
            self._headers = []

        def report(self, simulation, state):
            was_initialized = self._hasInitialized

            # spoof the file-like object with a string buffer
            self._out = StringIO()
            super(QueueStateDataReporter, self).report(simulation, state)

            if not was_initialized:
                # the first report has two lines on it -- we want to look at the first, as it contains the headers
                # print(self._out.getvalue())
                initial, line = self._out.getvalue().split('\n', 1)
                headers = initial.strip().split(',')
                # filter out some extra quotation marks and comment characters
                self._headers = [e.strip('#"\'') for e in headers]
            else:
                line = self._out.getvalue()
            t = [e.strip('%"\'') for e in line.strip().split(',')]
            # split the line based on whatever separator we know that the parent was using, and then cast to float

            msg = dict(zip(self._headers, t))
            queue.put(msg)

    return QueueStateDataReporter


##############################################################################
# Classes
##############################################################################

class Communicate(QtCore.QObject):
    dataSignal = QtCore.Signal(dict)
    main_signal = QtCore.Signal(dict)
    thread_id_keeper = QtCore.Signal(int)
    decomp_process = QtCore.Signal(list)
    finish_alert = QtCore.Signal(str)
    inform_about_situation = QtCore.Signal(str)


class OpenMMScriptRunner(QtCore.QObject):
    plots_created = bool
    openmm_script_code = str
    status = str
    pid_idents = []
    Signals = Communicate()
    plotdata = dict
    global _stop_running

    def __init__(self, script):
        super(OpenMMScriptRunner, self).__init__()
        # self.plotdata = dict
        self.plots_created = False
        self.decomp_data = []

        global _stop_running

        _stop_running = False

        self.openmm_script_code = script
        q = queue.Queue()

        self.t1 = threading.Thread(target=self.run_openmm_script, args=(self.openmm_script_code, q), daemon=True)
        self.t2 = threading.Thread(target=self.queue_consumer, args=(q,), daemon=True)

        self.t1.start()
        self.t2.start()

    def stop_threads(self):
        # global _stop_running
        # _stop_running = True
        """
        T = self.t1.is_running()
        print(T)
        self.t1.stop()
        """
        print("NO STOP FUNCTION YET")

    def run_openmm_script(self, code, queue):

        def fix_code():
            itoks = tokenize.generate_tokens(StringIO(code).readline)

            def run():
                for toktype, toktext, (srow, scol), (erow, ecol), line in itoks:
                    if toktext == 'StateDataReporter':
                        toktext = '__queue_reporter_factory(__queue)'

                    yield toktype, toktext, (srow, scol), (erow, ecol), line

            return tokenize.untokenize(run()) + '__queue.put(None)'

        try:
            code = fix_code()
        except tokenize.TokenError:
            raise ValueError('The script has a syntax error!')

        exec(code, {'__queue': queue, '__queue_reporter_factory': queue_reporter_factory})

    def queue_consumer(self, q):
        self.status = 'Running...'

        while True:
            try:
                if _stop_running:
                    break
                msg = q.get_nowait()
                if msg is None:
                    break
                self.update_plot(msg)

            except queue.Empty:
                time.sleep(0.1)

        self.status = 'Done'

    def create_plots(self, keys):
        self.plotdata = dict(zip(keys, [[]] * len(keys)))
        # figure out which key will be the x axis
        if 'Step' not in keys:
            raise ValueError('The reporter has not step information, so there is no x-axis to plot graphs!')

    def decomp_progress(self, data):
        self.decomp_data.append(data)

    def update_plot(self, msg):
        if not self.plots_created and type(msg) == dict:
            self.create_plots(msg.keys())
            self.plots_created = True

        if type(msg) == dict:
            for k, v in msg.items():
                current = self.plotdata.get(k)
                self.plotdata.update({k: np.concatenate((current, v), axis=None)})
            self.Signals.dataSignal.emit(self.plotdata)

        if type(msg) == list:
            self.Signals.decomp_process.emit(msg)

        if type(msg) == str:
            if msg == "Progress Finished Succesfully :)":
                self.Signals.finish_alert.emit(msg)
            else:
                self.Signals.inform_about_situation.emit(msg)


class Graphs(QWidget):
    pg.setConfigOption('background', None)
    pg.setConfigOption('foreground', (197, 198, 199))
    pg.setConfigOptions(antialias=True)
    global current_step_keeper

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.win = pg.GraphicsLayoutWidget(show=False, title="RMSD Plotting Region")

        self.real_time_as_minute = []
        self.real_speed = []
        self.current_step_keeper = None
        # self.win = pg.GraphicsWindow(show=False, title="Basic plotting examples")
        # setting style sheet to the plot window
        self.win.setStyleSheet("border : 2px solid green; padding: -5px; border-radius: 10px; """)
        self.win.setWindowTitle('Real Time Simulation Monitoring')

        self.temperature_graph = self.win.addPlot(title="Temperature")
        self.temperature_graph.disableAutoRange()
        self.temperature_graph.addLegend()
        self.temperature_graph.getViewBox().setBackgroundColor((129, 105, 161, 20))
        self.temperature_graph.setLabel('left', "Temperature", units='K')
        self.temperature_graph.setLabel('bottom', "Step")
        # self.temperature_graph.setLogMode(x=True, y=False) #logaritmik mode
        self.temperature_graph.setLogMode(x=False, y=False)
        self.temperature_graph.setYRange(200, 400, padding=0)
        self.temperature_graph.showGrid(x=True, y=True)
        self.temperature_graph.getAxis('left').enableAutoSIPrefix(False)
        self.temperature_graph_plot = self.temperature_graph.plot(name='Temperature')

        # self.win.nextRow()
        self.energy_graph = self.win.addPlot(title="Energy")
        self.energy_graph.addLegend()
        self.energy_graph.getViewBox().setBackgroundColor((129, 105, 161, 20))
        self.energy_graph.setLabel('left', "Energy", units='kJ/mole')
        self.energy_graph.setLabel('bottom', "Step")
        # self.temperature_graph.setLogMode(x=True, y=False) #logaritmik mode
        self.energy_graph.showGrid(x=True, y=True)

        self.potential_energy_graph = self.energy_graph.plot(name='Potential')
        self.kinetic_energy_graph = self.energy_graph.plot(name='Kinetic')
        self.total_energy_graph = self.energy_graph.plot(name='Total')
        # self.energy_graph.addLegend()

        self.win.nextRow()
        self.simulation_speed_graph = self.win.addPlot(title="Speed", row=1, colspan=2)
        self.simulation_speed_graph.addLegend()
        self.simulation_speed_graph.getViewBox().setBackgroundColor((129, 105, 161, 20))
        self.simulation_speed_graph.setLabel('left', "Speed", units='ns/day')
        self.simulation_speed_graph.setLabel('bottom', "Step")
        self.simulation_speed_graph.showGrid(x=True, y=True)

        self.simulation_speed_graph_plot = self.simulation_speed_graph.plot(name='Speed (ns/day)')

        self.win.nextRow()
        self.simulation_time_graph = self.win.addPlot(title="Remaining Time", row=2, colspan=2)
        self.simulation_time_graph.addLegend()
        self.simulation_time_graph.getViewBox().setBackgroundColor((129, 105, 161, 20))
        self.simulation_time_graph.setLabel('left', "Remaining Time", units='sec')
        self.simulation_time_graph.setLabel('bottom', "Step")
        self.simulation_time_graph.showGrid(x=True, y=True)

        self.simulation_time_graph_plot = self.simulation_time_graph.plot(name='Remaining Time (sec)')
        self.simulation_time_graph.setLogMode(x=False, y=False)  # logaritmik mode
        self.simulation_time_graph.enableAutoRange(axis='y')

        self.first_entrance = 0

    def pretty_speed(self, ins_speed):
        """Format the speed (ns/day) as pretty"""
        speed_style = ins_speed.split(':')
        if ins_speed == '':
            return self.real_speed.append(float(0))

        if len(speed_style) == 1:

            if speed_style[0] == '--':
                return self.real_speed.append(float(0))
            return self.real_speed.append(float(speed_style[0]))

        if len(speed_style) == 2:
            if speed_style[0] == '--' or speed_style[1] == '--':
                return self.real_speed.append(float(0))
            return self.real_speed.append(float(ins_speed))

    def pretty_time(self, t_remaining):
        """Format the time as minute"""
        time_style = t_remaining.split(':')
        print(t_remaining)
        if t_remaining == '':
            return self.real_time_as_minute.append(float(0))

        if len(time_style) == 1:
            if time_style[0] == '--':
                return self.real_time_as_minute.append(float(0))
            second = float(time_style[1]) / 60
            return self.real_time_as_minute.append(second)

        if len(time_style) == 2:
            if time_style[0] == '--' or time_style[1] == '--':
                return self.real_time_as_minute.append(float(0))
            minute = float(time_style[0])
            second = float(time_style[1]) / 60
            return self.real_time_as_minute.append(minute + second)

        if len(time_style) == 3:
            if time_style[0] == '--' or time_style[1] == '--' or time_style[2] == '--':
                return self.real_time_as_minute.append(float(0))
            hour = float(time_style[0]) * 60
            minute = float(time_style[1])
            second = float(time_style[2]) / 60
            return self.real_time_as_minute.append(hour + minute + second)

    def update_graph(self, data):
        x = np.array(data["Step"], dtype=np.float)
        y_temp = np.array(data["Temperature (K)"], dtype=np.float)
        y_potential = np.array(data["Potential Energy (kJ/mole)"], dtype=np.float)
        y_kinetic = np.array(data["Kinetic Energy (kJ/mole)"], dtype=np.float)
        y_total = np.array(data["Total Energy (kJ/mole)"], dtype=np.float)

        y_speed = np.array(data["Speed (ns/day)"])[-1]
        y_time_remaining = np.array(data["Time Remaining"])[-1]

        self.pretty_speed(y_speed)
        self.pretty_time(y_time_remaining)

        if x.shape == y_temp.shape:

            if self.current_step_keeper is not None and self.current_step_keeper[-1] > x[-1]:
                if self.current_step_keeper[-1] - self.current_step_keeper[-2] > 1:

                    ticks = pg.VTickGroup(xvals=[self.current_step_keeper[-1]], yrange=[0, 2.5],
                                          pen={'color': 'g', 'width': 2.5})
                    self.simulation_time_graph_plot.getViewBox().addItem(ticks)
                x = np.append(self.current_step_keeper, self.current_step_keeper[-1] + 1)

            try:
                self.temperature_graph_plot.setData(x=x, y=y_temp, clear=True, pen=pg.mkPen((255, 0, 0), width=3),
                                                    name="Temperature", fillLevel=0.0, brush=(150, 150, 50, 30))
                # self.temperature_graph.autoRange()

                self.potential_energy_graph.setData(x=x, y=y_potential, clear=True, pen=pg.mkPen((255, 0, 0), width=3),
                                                    name="Potential")

                self.kinetic_energy_graph.setData(x=x, y=y_kinetic, clear=True, pen=pg.mkPen((0, 255, 0), width=3),
                                                  name="Kinetic")

                self.total_energy_graph.setData(x=x, y=y_total, clear=True, pen=pg.mkPen((0, 0, 255), width=3),
                                                fillLevel=0.0, brush=(150, 150, 50, 10), name="Total")

                self.simulation_time_graph_plot.setData(x=x, y=self.real_time_as_minute,
                                                        pen=pg.mkPen((0, 0, 255), width=3),
                                                        fillLevel=0.0, name="Rime Remaining (sec)",
                                                        brush=(150, 150, 50, 10))

                self.simulation_speed_graph_plot.setData(x=x, y=self.real_speed, pen=pg.mkPen((200, 200, 200), width=3),
                                                         symbolBrush=(255, 0, 0), symbolPen='w', fillLevel=0.0,
                                                         name="Speed",
                                                         brush=(150, 150, 50, 30))

                self.current_step_keeper = x

            except Exception as err:
                print("========================\n", err)
                print(x, "\n", self.real_speed)
                print(y_speed)
                print("\n========================")

    def updating_decomposion(self, data_decomp):
        print(data_decomp)

    def run_script(self, contents):
        self.contents = contents
        self.runner = OpenMMScriptRunner(self.contents)
        self.runner.Signals.dataSignal.connect(lambda plotdata: self.update_graph(plotdata))
        # self.runner.Signals.decomp_process.connect(lambda decomp_data: self.updating_decomposion(decomp_data))

    def stop_th(self):
        self.runner.stop_threads()

# if __name__ == '__main__':
#     app = QtWidgets.QApplication(sys.argv)
#     main = Graphs(contents)
#     main.show()
#     sys.exit(app.exec_())
