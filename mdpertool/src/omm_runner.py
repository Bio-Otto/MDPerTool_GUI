import os
import re
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
import subprocess


# ################################################################################################################### #
# #################################################### FUNCTIONS #################################################### #
# ################################################################################################################### #

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
    run_speed = QtCore.Signal(float)
    finish_alert = QtCore.Signal(str)
    inform_about_situation = QtCore.Signal(str)


class OpenMMScriptRunner(QtCore.QObject):
    plots_created = bool
    openmm_script_code = str
    status = str
    pid_idents = []
    Signals = Communicate()
    plotdata = dict
    process = None
    decomp_data = []
    speed_data = []
    def __init__(self, script):
        super(OpenMMScriptRunner, self).__init__()
        # self.plotdata = dict
        self.plots_created = False
        #self.decomp_data = []
        #self.speed_data = []
        self.process = None

        self.openmm_script_code = script
        q = queue.Queue()

        self.t1 = threading.Thread(target=self.run_openmm_script, args=(self.openmm_script_code, q), daemon=True)
        self.t2 = threading.Thread(target=self.queue_consumer, args=(q,), daemon=True)

        self.t1.start()
        self.t2.start()

    def stop_threads(self):
        self.stop_process()
        self.plots_created = False

    def run_openmm_script(self, code, queue):
        def fix_code():
            itoks = tokenize.generate_tokens(StringIO(code).readline)

            def run():
                for toktype, toktext, (srow, scol), (erow, ecol), line in itoks:
                    yield toktype, toktext, (srow, scol), (erow, ecol), line

            return tokenize.untokenize(run())

        try:
            code = fix_code()
        except tokenize.TokenError:
            raise ValueError('The script has a syntax error!')

        with open('temp_script.py', 'w') as f:
            f.write(code)

        self.process = subprocess.Popen(['python', 'temp_script.py'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        while True:
            output = self.process.stdout.readline()
            if output == b'' and self.process.poll() is not None:
                break
            if output:
                queue.put(output.decode().strip())

        while True:
            error_output = self.process.stderr.readline()
            if error_output == b'' and self.process.poll() is not None:
                break
            if error_output:
                queue.put(error_output.decode().strip())

        self.process.wait()

        if self.process.returncode != 0:
            raise ValueError('Script execution failed!')

        os.remove('temp_script.py')

    # =============================================

    def stop_process(self):
        try:
            self.process.terminate()
            self.process.wait()
        except Exception as err:
            print("ERROR IN STOP PROCESS: %s" % err)

    def queue_consumer(self, q):
        self.status = 'Running...'
        _headers = []
        decompose_started = None

        while True:
            try:

                msg = q.get_nowait()
                if msg is None:
                    break

                else:
                    pass
                    # print("======================================================================")
                    # print("MESSAGE: %s" % msg)
                    # print("TYPE: %s" % type(msg))
                    # print("======================================================================")

                if 'INFO |' in msg:
                    info_log = re.search(r'INFO \| (.+)', msg)
                    self.update_plot(info_log.group(1))

                elif 'CRITICAL |' in msg:
                    critic_log = re.search(r'CRITICAL \| (.+)', msg)
                    self.update_plot(critic_log.group(1))

                elif 'WARNING |' in msg:
                    warning_log = re.search(r'WARNING \| (.+)', msg)
                    self.update_plot(warning_log.group(1))

                elif 'ERROR |' in msg:
                    error_log = re.search(r'ERROR \| (.+)', msg)
                    self.update_plot(error_log.group(1))

                elif '#"Progress (%)"' in msg:
                    # the first report has two lines on it -- we want to look at the first, as it contains the headers
                    # print(self._out.getvalue())
                    headers = msg.strip().split(',')
                    # filter out some extra quotation marks and comment characters
                    _headers = [e.strip('#"\'') for e in headers]

                elif type(msg) is not dict:
                    t = [e.strip('%"\'') for e in msg.strip().split(',')]
                    if len(_headers) == len(t):
                        msg = dict(zip(_headers, t))
                        q.put(msg)
                        self.update_plot(msg)

                if 'DECOMPOSE |' in msg:
                    decomp_info_log = re.search(r"Decomposition Progress: ([\d.]+)", msg)

                    if decomp_info_log:
                        decompose_started = True
                        extracted_number = float(decomp_info_log.group(1))
                        formatted_number = round(extracted_number, 2)
                        self.decomp_data.append(formatted_number)
                        #self.Signals.decomp_process.emit(self.decomp_data)
                        self.update_plot(self.decomp_data)


                """
                elif type(msg) is not dict:
                    if decompose_started:
                        print("NEW MESSAGE: %s" % msg)
                        break

                    if 'INFO | Decompose started using XTC File' in msg:
                        print("passed to decompose.")
                        decompose_started = True
                        break

                    t = [e.strip('%"\'') for e in msg.strip().split(',')]
                    # split the line based on whatever separator we know that the parent was using, and then cast to
                    # float
                    if len(_headers) > 0:
                        msg = dict(zip(_headers, t))
                        q.put(msg)

                    print(msg)
                    self.update_plot(msg)
                """
            except queue.Empty:
                time.sleep(0.1)

        self.status = 'Done'

    def create_plots(self, keys):
        self.plotdata = dict(zip(keys, [[]] * len(keys)))
        # figure out which key will be the x-axis
        if 'Step' not in keys:
            raise ValueError('The reporter has not step information, so there is no x-axis to plot graphs!')

    def speed_reporter(self, data):
        self.speed_data.append(data)

    def update_plot(self, msg):

        if self.plots_created and type(msg) == dict:
            for k, v in msg.items():
                current = self.plotdata.get(k)

                self.plotdata.update({k: np.concatenate((current, v), axis=None)})
            self.Signals.dataSignal.emit(self.plotdata)

        if not self.plots_created and type(msg) == dict:
            self.create_plots(msg.keys())
            self.plots_created = True

        if type(msg) == list:
            self.Signals.decomp_process.emit(msg)

        if type(msg) == float:
            print("here speed")
            self.Signals.run_speed.emit(msg)

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
        self.win.setStyleSheet("border : 2px solid green; padding: -5px; border-radius: 10px; """)
        self.win.setWindowTitle('Real Time Simulation Monitoring')

        self.temperature_graph = self.win.addPlot(title="Temperature")

        # self.temperature_graph.disableAutoRange()
        self.temperature_graph.addLegend()
        self.temperature_graph.getViewBox().setBackgroundColor((129, 105, 161, 20))
        self.temperature_graph.setLabel('left', "Temperature", units='K')
        self.temperature_graph.setLabel('bottom', "Step")
        # self.temperature_graph.setLogMode(x=True, y=False) #logaritmik mode
        self.temperature_graph.setLogMode(x=False, y=False)
        # self.temperature_graph.setYRange(200, 400, padding=0)
        self.temperature_graph.showGrid(x=True, y=True)
        self.temperature_graph.getAxis('left').enableAutoSIPrefix(False)
        self.temperature_graph_plot = pg.PlotDataItem(clear=True, pen=pg.mkPen((255, 0, 0), width=3),
                                                      name="Temperature", fillLevel=0.0, brush=(150, 150, 50, 30))
        self.temperature_graph.addItem(self.temperature_graph_plot)
        # self.temperature_graph_plot = self.temperature_graph.plot(name='Temperature')

        # self.win.nextRow()
        self.energy_graph = self.win.addPlot(title="Energy")
        # self.energy_graph.disableAutoRange()
        self.energy_graph.addLegend()
        self.energy_graph.getViewBox().setBackgroundColor((129, 105, 161, 20))
        self.energy_graph.setLabel('left', "Energy", units='kJ/mole')
        self.energy_graph.setLabel('bottom', "Step")
        # self.temperature_graph.setLogMode(x=True, y=False) #logaritmik mode
        self.energy_graph.setLogMode(x=False, y=False)
        # self.energy_graph.setYRange(200, 400, padding=0)
        self.energy_graph.showGrid(x=True, y=True)
        self.energy_graph.getAxis('left').enableAutoSIPrefix(False)

        self.potential_energy_graph = pg.PlotDataItem(clear=True, pen=pg.mkPen((255, 0, 0), width=2),
                                                      name="Potential", fillLevel=0.0, brush=(150, 150, 50, 30))
        self.kinetic_energy_graph = pg.PlotDataItem(clear=True, pen=pg.mkPen((0, 255, 0), width=2),
                                                    name="Kinetic", fillLevel=0.0, brush=(150, 150, 50, 30))
        self.total_energy_graph = pg.PlotDataItem(clear=True, pen=pg.mkPen((0, 0, 255), width=2),
                                                  name="Total", fillLevel=0.0, brush=(150, 150, 50, 30))
        self.energy_graph.addItem(self.potential_energy_graph)
        self.energy_graph.addItem(self.kinetic_energy_graph)
        self.energy_graph.addItem(self.total_energy_graph)
        # self.energy_graph.addLegend()

        self.win.nextRow()
        self.simulation_speed_graph = self.win.addPlot(title="Speed", row=1, colspan=2)
        # self.simulation_speed_graph.disableAutoRange()
        self.simulation_speed_graph.addLegend()
        self.simulation_speed_graph.getViewBox().setBackgroundColor((129, 105, 161, 20))
        self.simulation_speed_graph.setLabel('left', "Speed", units='ns/day')
        self.simulation_speed_graph.setLabel('bottom', "Step")
        self.simulation_speed_graph.setLogMode(x=False, y=False)
        # self.simulation_speed_graph.setYRange(200, 400, padding=0)
        self.simulation_speed_graph.showGrid(x=True, y=True)
        self.simulation_speed_graph.getAxis('left').enableAutoSIPrefix(False)

        self.simulation_speed_graph_plot = pg.PlotDataItem(clear=True, name="Speed (ns/day)", fillLevel=0.0,
                                                           brush=(150, 150, 50, 30),
                                                           pen=pg.mkPen((200, 200, 200), width=2),
                                                           symbolBrush=(255, 0, 0), symbolPen='w',
                                                           )
        self.simulation_speed_graph.addItem(self.simulation_speed_graph_plot)

        self.win.nextRow()
        self.simulation_time_graph = self.win.addPlot(title="Remaining Time", row=2, colspan=2)
        # self.simulation_time_graph.disableAutoRange()
        self.simulation_time_graph.addLegend()
        self.simulation_time_graph.getViewBox().setBackgroundColor((129, 105, 161, 20))
        self.simulation_time_graph.setLabel('left', "Remaining Time", units='sec')
        self.simulation_time_graph.setLabel('bottom', "Step")
        self.simulation_time_graph.setLogMode(x=False, y=False)
        # self.simulation_time_graph.setYRange(200, 400, padding=0)
        self.simulation_time_graph.showGrid(x=True, y=True)
        # self.simulation_time_graph.setLogMode(x=False, y=False)  # logaritmik mode
        # self.simulation_time_graph.enableAutoRange(axis='y')
        self.simulation_time_graph.getAxis('left').enableAutoSIPrefix(False)

        self.simulation_time_graph_plot = pg.PlotDataItem(clear=True, pen=pg.mkPen((255, 255, 0), width=2),
                                                          fillLevel=0.0, name="Remaining Time (sec)",
                                                          brush=(150, 150, 50, 10))
        self.simulation_time_graph.addItem(self.simulation_time_graph_plot)
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
            second = float(time_style[0]) / 60
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
        print("TIME REMAINING: %s" % t_remaining)

    def update_graph(self, data):

        try:
            x = np.array(data["Step"], dtype=np.float64)
            y_temp = np.array(data["Temperature (K)"], dtype=np.float64)
            y_potential = np.array(data["Potential Energy (kJ/mole)"], dtype=np.float64)
            y_kinetic = np.array(data["Kinetic Energy (kJ/mole)"], dtype=np.float64)
            y_total = np.array(data["Total Energy (kJ/mole)"], dtype=np.float64)

            y_speed = np.array(data["Speed (ns/day)"])[-1]
            self.pretty_speed(y_speed)

            y_time_remaining = np.array(data["Time Remaining"])[-1]
            self.pretty_time(y_time_remaining)

            if x.shape == y_temp.shape:

                if self.current_step_keeper is not None and self.current_step_keeper[-1] > x[-1]:
                    if self.current_step_keeper[-1] - self.current_step_keeper[-2] > 1:
                        ticks = pg.VTickGroup(xvals=[self.current_step_keeper[-1]], yrange=[0, 2.5],
                                              pen={'color': 'g', 'width': 2.5})
                        self.simulation_time_graph_plot.getViewBox().addItem(ticks)
                    x = np.append(self.current_step_keeper, self.current_step_keeper[-1] + 1)

                try:

                    self.temperature_graph_plot.setData(x=x, y=y_temp)
                    # self.temperature_graph.autoRange()
                    self.potential_energy_graph.setData(x=x, y=y_potential)
                    self.kinetic_energy_graph.setData(x=x, y=y_kinetic)
                    self.total_energy_graph.setData(x=x, y=y_total)

                    # self.simulation_time_graph_plot.setData(x=x, y=self.real_time_as_minute)
                    # self.simulation_speed_graph_plot.setData(x=x, y=self.real_speed)

                    self.current_step_keeper = x

                except Exception as err:
                    import traceback
                    print("ERROR 1: %s" % err)
                    traceback.print_exc()

        except Exception as e:
            print("ERROR 2: %s" % e)

    def updating_decomposion(self, data_decomp):
        pass

    def updating_current_speed(self, data_speed):
        print("SPEEEEEEDDDDD: ", data_speed)
        pass

    def run_script(self, contents):
        self.contents = contents
        self.runner = OpenMMScriptRunner(self.contents)
        self.runner.Signals.dataSignal.connect(lambda plotdata: self.update_graph(plotdata))
        #self.runner.Signals.decomp_process.connect(lambda decomp_data: self.updating_decomposion(decomp_data))
        self.runner.Signals.run_speed.connect(lambda speed_data: self.updating_current_speed(speed_data))

    def stop_th(self):
        self.runner.stop_threads()

# if __name__ == '__main__':
#     app = QtWidgets.QApplication(sys.argv)
#     main = Graphs(contents)
#     main.show()
#     sys.exit(app.exec_())
