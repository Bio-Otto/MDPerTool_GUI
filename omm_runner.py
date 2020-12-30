from simtk.openmm.app import StateDataReporter
from io import StringIO
import time
import ctypes
import queue
import threading
import multiprocessing
import tokenize
import pystache
import numpy as np
from PyQt5 import QtWidgets, uic, QtCore, QtGui
from pyqtgraph import PlotWidget, plot, dockarea
import pyqtgraph as pg
from PyQt5.QtCore import QTimer, QDateTime, pyqtSlot
import subprocess
import continuous_threading
import os
import sys
from PyQt5 import QtCore
from PyQt5.QtCore import (QCoreApplication, QPropertyAnimation, QDate, QDateTime, QMetaObject, QObject, QPoint, QRect,
                          QSize, QTime, QUrl, Qt, QEvent)
from PyQt5.QtGui import (QBrush, QColor, QConicalGradient, QCursor, QFont, QFontDatabase, QIcon, QKeySequence,
                         QLinearGradient, QPalette, QPainter, QPixmap, QRadialGradient)
from PyQt5.QtWidgets import *


# ##############################################################################
# # Functions
# ##############################################################################

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
                print(self._out.getvalue())
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
    dataSignal = QtCore.pyqtSignal(dict)
    main_signal = QtCore.pyqtSignal(dict)
    thread_id_keeper = QtCore.pyqtSignal(int)


class OpenMMScriptRunner(QtCore.QObject):
    plots_created = bool
    openmm_script_code = str
    status = str
    pid_idents = []
    Signals = Communicate()
    global _stop_running

    def __init__(self, script):
        super(OpenMMScriptRunner, self).__init__()
        self.plotdata = dict
        self.plots_created = False

        global _stop_running

        _stop_running = False

        self.openmm_script_code = script
        q = queue.Queue()

        self.t1 = threading.Thread(target=self.run_openmm_script, args=(self.openmm_script_code, q), daemon=True)
        self.t2 = threading.Thread(target=self.queue_consumer, args=(q,), daemon=True)

        self.t1.start()
        self.t2.start()

        #

    def stop_threads(self):

        # global _stop_running
        # _stop_running = True
        """
        T = self.t1.is_running()
        print(T)
        self.t1.stop()
        """

        ctypes.pythonapi.PyThreadState_SetAsyncExc(threading.get_ident(), ctypes.py_object(SystemExit))


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

    def update_plot(self, msg):
        if not self.plots_created:
            self.create_plots(msg.keys())
            self.plots_created = True

        for k, v in msg.items():
            current = self.plotdata.get(k)
            self.plotdata.update({k: np.concatenate((current, v), axis=None)})
        self.Signals.dataSignal.emit(self.plotdata)


# ###########################################################################

class Graphs(QWidget):
    # global curve, data, p6
    def __init__(self, *args, **kwargs):
        # self.contents = contents
        # super().__init__()
        super().__init__(*args, **kwargs)

        self.real_time_as_minute = []
        self.real_speed = []
        pg.setConfigOption('background', None)
        pg.setConfigOption('foreground', (197, 198, 199))
        self.win = pg.GraphicsWindow(show=True, title="Basic plotting examples")
        # setting style sheet to the plot window
        self.win.setStyleSheet("border : 2px solid green; padding: -5px; border-radius: 10px; """)
        self.win.setWindowTitle('Real Time Simulation Monitoring')
        # self.verticalLayout_16.addWidget(self.win)

        self.temperature_graph = self.win.addPlot(title="Temperature")
        self.temperature_graph.addLegend()
        # self.temperature_graph.showLabel('bottom', show=True)
        self.temperature_graph.getViewBox().setBackgroundColor((129, 105, 161, 20))
        self.temperature_graph.setLabel('left', "Temperature", units='K')
        self.temperature_graph.setLabel('bottom', "Step")
        # self.temperature_graph.setLogMode(x=True, y=False) #logaritmik mode
        self.temperature_graph.setYRange(200, 400, padding=0)
        self.temperature_graph.showGrid(x=True, y=True)

        # self.win.nextRow()
        self.energy_graph = self.win.addPlot(title="Energy")

        self.energy_graph.getViewBox().setBackgroundColor((129, 105, 161, 20))
        self.energy_graph.setLabel('left', "Energy", units='kJ/mole')
        self.energy_graph.setLabel('bottom', "Step")
        # self.temperature_graph.setLogMode(x=True, y=False) #logaritmik mode
        self.energy_graph.showGrid(x=True, y=True)

        self.potential_energy_graph = self.energy_graph.plot()
        self.kinetic_energy_graph = self.energy_graph.plot()
        self.total_energy_graph = self.energy_graph.plot()
        self.energy_graph.addLegend()

        self.win.nextRow()
        self.simulation_speed_and_time_graph = self.win.addPlot(title="Speed")
        self.simulation_speed_and_time_graph.addLegend()

        self.simulation_speed_graph = self.simulation_speed_and_time_graph.plot()
        self.simulation_time_graph = self.simulation_speed_and_time_graph.plot()

        self.simulation_speed_and_time_graph.getViewBox().setBackgroundColor((129, 105, 161, 20))
        self.simulation_speed_and_time_graph.setLabel('left', "Speed", units='ns/day')
        self.simulation_speed_and_time_graph.setLabel('bottom', "Step")
        self.simulation_speed_and_time_graph.showGrid(x=True, y=True)

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
            self.temperature_graph.plot(x=x, y=y_temp, clear=True, pen=pg.mkPen((255, 0, 0), width=3),
                                        name="Temperature", fillLevel=0.0, brush=(150, 150, 50, 30))
            self.temperature_graph.autoRange()

            self.potential_energy_graph.setData(x=x, y=y_potential, clear=True, pen=pg.mkPen((255, 0, 0), width=3),
                                                name="Potential")

            self.kinetic_energy_graph.setData(x=x, y=y_kinetic, clear=True, pen=pg.mkPen((0, 255, 0), width=3),
                                              name="Kinetic")

            self.total_energy_graph.setData(x=x, y=y_total, clear=True, pen=pg.mkPen((0, 0, 255), width=3),
                                            fillLevel=0.0, brush=(150, 150, 50, 10), name="Total")

            self.simulation_time_graph.setData(x=x, y=self.real_time_as_minute, pen=pg.mkPen((0, 0, 255), width=3),
                                               fillLevel=0.0, name="Rime Remaining (sec)", brush=(150, 150, 50, 10))

            self.simulation_speed_graph.setData(x=x, y=self.real_speed, pen=pg.mkPen((200, 200, 200), width=3),
                                                symbolBrush=(255, 0, 0), symbolPen='w', fillLevel=0.0, name="Speed",
                                                brush=(150, 150, 50, 30))

            # if self.first_entrance == 1:
            #     ay = self.simulation_speed_and_time_graph.getAxis('left')
            #     dy = [(value, '{:.3f}'.format(value)) for value in y_speed]
            #     ay.setTicks([dy, []])
            #     print(y_speed)
            # QtCore.QCoreApplication.processEvents()
            # p2.plot(np.random.normal(size=120) + 10, pen=(0, 0, 255), name="Blue curve")

    def run_script(self, contents):
        self.contents = contents
        self.runner = OpenMMScriptRunner(self.contents)
        self.runner.Signals.dataSignal.connect(lambda plotdata: self.update_graph(plotdata))

    def stop_th(self):
        self.runner.stop_threads()

# if __name__ == '__main__':
#     app = QtWidgets.QApplication(sys.argv)
#     main = Graphs(contents)
#     main.show()
#     sys.exit(app.exec_())
