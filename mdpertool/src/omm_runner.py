import os
import re
import queue
import threading
import time
import tokenize
import traceback
import logging
import sys
from io import StringIO
import numpy as np
import pyqtgraph as pg
from PySide2 import QtCore
from PySide2.QtWidgets import QMessageBox, QWidget
from openmm.app import StateDataReporter
import subprocess
import tempfile
from datetime import datetime
from .message import Message_Boxes
from gui.ui_styles import Style


##############################################################################
# Classes
##############################################################################

class Communicate(QtCore.QObject):
    dataSignal = QtCore.Signal(dict)
    main_signal = QtCore.Signal(dict)
    thread_id_keeper = QtCore.Signal(int)
    decomp_process = QtCore.Signal(list)
    classic_md_remain_time = QtCore.Signal(list)
    reference_md_remain_time = QtCore.Signal(list)
    dissipation_md_remain_time = QtCore.Signal(list)
    run_speed = QtCore.Signal(list)
    finish_alert = QtCore.Signal(str)
    inform_about_situation = QtCore.Signal(str)
    real_time_pertub_info = QtCore.Signal(str)


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
    classic_md_remain_data = []
    reference_md_remain_data = []
    dissipation_md_remain_data = []
    pymol_statu = str


    def __init__(self, script):
        super(OpenMMScriptRunner, self).__init__()
        self.plots_created = False
        self.process = None
        self.stop_requested = False
        self.shutdown_requested = threading.Event()

        self.openmm_script_code = script
        self._queue = queue.Queue()
        self.max_plot_points = 3000
        self.last_plot_emit_time = 0.0
        self.plot_emit_interval_sec = 0.15

        self.t1 = threading.Thread(target=self.run_openmm_script, args=(self.openmm_script_code, self._queue),
                                   daemon=True)
        self.t2 = threading.Thread(target=self.queue_consumer, args=(self._queue,), daemon=True)

        self.t1.start()
        self.t2.start()

    def _safe_emit(self, signal_obj, payload):
        if self.shutdown_requested.is_set():
            return False
        try:
            signal_obj.emit(payload)
            return True
        except RuntimeError:
            self.shutdown_requested.set()
            return False

    def stop_threads(self):
        self.shutdown_requested.set()
        was_stopped = self.stop_process()
        if hasattr(self, '_queue'):
            try:
                self._queue.put(None)
            except Exception:
                pass
        if hasattr(self, 't2') and self.t2.is_alive():
            self.t2.join(timeout=2)
        if was_stopped:
            self.plots_created = False
        return was_stopped

    def shutdown(self):
        self.shutdown_requested.set()
        try:
            if self.process is not None and self.process.poll() is None:
                self.stop_requested = True
                self.process.terminate()
                self.process.wait(timeout=3)
        except Exception:
            try:
                if self.process is not None and self.process.poll() is None:
                    self.process.kill()
            except Exception:
                pass

        if hasattr(self, '_queue'):
            try:
                self._queue.put(None)
            except Exception:
                pass

        if hasattr(self, 't2') and self.t2.is_alive():
            self.t2.join(timeout=2)

    def run_openmm_script(self, code, queue):
        def _safe_decode(raw_line):
            if isinstance(raw_line, str):
                return raw_line
            for encoding in ('utf-8', 'cp1254', 'latin-1'):
                try:
                    return raw_line.decode(encoding)
                except UnicodeDecodeError:
                    continue
            return raw_line.decode('utf-8', errors='replace')

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

        with tempfile.NamedTemporaryFile(mode='w', suffix='.py', delete=False, encoding='utf-8') as temp_script:
            temp_script.write(code)
            temp_script_path = temp_script.name

        debug_script_path = None
        try:
            debug_dir = os.path.join(os.getcwd(), 'output', 'debug_scripts')
            os.makedirs(debug_dir, exist_ok=True)
            stamp = datetime.now().strftime('%Y%m%d_%H%M%S')
            debug_script_path = os.path.join(debug_dir, f'generated_openmm_script_{stamp}.py')
            with open(debug_script_path, 'w', encoding='utf-8') as debug_script_file:
                debug_script_file.write(code)
            queue.put(f"INFO | Debug script saved: {debug_script_path}")
        except Exception as debug_file_error:
            queue.put(f"WARNING | Unable to save debug script file: {debug_file_error}")

        try:
            child_env = os.environ.copy()
            child_env['PYTHONIOENCODING'] = 'utf-8'
            child_env['PYTHONUTF8'] = '1'

            self.process = subprocess.Popen(
                [sys.executable, '-u', temp_script_path],
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                bufsize=0,
                env=child_env
            )

            for output in iter(self.process.stdout.readline, b''):
                if output:
                    queue.put(_safe_decode(output).rstrip('\r\n'))

            self.process.wait()

            if self.process.returncode != 0:
                if self.stop_requested:
                    queue.put("INFO | Run stopped by user.")
                else:
                    queue.put(f"ERROR | Script process exited with code {self.process.returncode}.")
                    queue.put("Script execution failed. Please review the traceback/messages above.")
                    if debug_script_path:
                        queue.put(f"ERROR | Re-run this file directly for debugging: {debug_script_path}")
                # raise subprocess.CalledProcessError(self.process.returncode, 'Script execution has been stoped by the '
                #                                                             'User!')

        except (subprocess.CalledProcessError, OSError, Exception) as e:
            logging.exception("Exception while running generated OpenMM script")
            queue.put(f"CRITICAL | Runner exception: {e}")
            queue.put(traceback.format_exc())

        finally:
            self.stop_requested = False
            if os.path.exists(temp_script_path):
                os.remove(temp_script_path)
            queue.put(None)

    def stop_process(self):

        if self.process is not None and self.process.poll() is None:
            try:
                """Close Application Question Message Box."""
                close_answer = Message_Boxes.Question_message(self, "Are you sure!", "Do you really want to stop the "
                                                                                     "run?",
                                                              Style.MessageBox_stylesheet)

                if close_answer == QMessageBox.Yes:
                    try:
                        self.stop_requested = True
                        self.process.terminate()
                        self.process.wait()
                        return True
                    except (subprocess.CalledProcessError, OSError) as err:
                        print("ERROR IN STOP PROCESS: %s" % err)
                        return False

                if close_answer == QMessageBox.No:
                    return False

            except Exception as inst:
                Message_Boxes.Critical_message(self, 'An unexpected error has occurred!', str(inst),
                                               Style.MessageBox_stylesheet)
                return False
        else:
            Message_Boxes.Information_message(self, "Info", "There is no an active run", Style.MessageBox_stylesheet)
            return False

    def queue_consumer(self, q):
        self.status = 'Running...'
        _headers = []
        decompose_started = None

        while not self.shutdown_requested.is_set():
            try:

                msg = q.get_nowait()
                if msg is None:
                    break

                msg_handled = False

                if type(msg) == str and 'INFO |' in msg:
                    info_log = re.search(r'INFO \|\s*(.*)', msg)
                    parsed_info = info_log.group(1) if info_log else msg
                    self.update_plot(parsed_info)
                    msg_handled = True

                elif type(msg) == str and 'CRITICAL |' in msg:
                    critic_log = re.search(r'CRITICAL \|\s*(.*)', msg)
                    parsed_critical = critic_log.group(1) if critic_log else msg
                    self.update_plot(parsed_critical)
                    msg_handled = True

                elif type(msg) == str and 'WARNING |' in msg:
                    warning_log = re.search(r'WARNING \|\s*(.*)', msg)
                    parsed_warning = warning_log.group(1) if warning_log else msg
                    self.update_plot(parsed_warning)
                    msg_handled = True

                elif type(msg) == str and 'ERROR |' in msg:
                    error_log = re.search(r'ERROR \|\s*(.*)', msg)
                    parsed_error = error_log.group(1) if error_log else msg
                    print(parsed_error)
                    self.update_plot(parsed_error)
                    msg_handled = True

                elif type(msg) == str and 'Progress (%)' in msg and '|' not in msg:
                    headers = msg.strip().split(',')
                    _headers = [e.strip('#"\'') for e in headers]
                    if 'Dissipation MD Progress (%)' in msg:
                        self.pymol_statu = "Started on PyMOL"
                        self._safe_emit(self.Signals.real_time_pertub_info, self.pymol_statu)
                    msg_handled = True

                elif type(msg) is not dict:
                    t = [e.strip('%"\'') for e in msg.strip().split(',')]

                    if len(_headers) == len(t):
                        msg = dict(zip(_headers, t))
                        q.put(msg)
                        self.update_plot(msg)
                        msg_handled = True

                if type(msg) == str and 'DECOMPOSE |' in msg:
                    decomp_info_log = re.search(r"Decomposition Progress: ([\d.]+)", msg)

                    if decomp_info_log:
                        decompose_started = True
                        extracted_number = float(decomp_info_log.group(1))
                        formatted_number = round(extracted_number, 2)
                        self.decomp_data.append(formatted_number)
                        self.update_plot(self.decomp_data)

                    self.pymol_statu = "Finished on PyMOL"
                    self._safe_emit(self.Signals.real_time_pertub_info, self.pymol_statu)
                    msg_handled = True

                if type(msg) == str and 'AUTO_EXTEND |' in msg:
                    self.pymol_statu = "Auto-Extending"
                    self._safe_emit(self.Signals.real_time_pertub_info, self.pymol_statu)
                    msg_handled = True

                if type(msg) == str and not msg_handled:
                    stripped_msg = msg.strip()
                    if stripped_msg:
                        self.update_plot(stripped_msg)

            except queue.Empty:
                time.sleep(0.1)
            except RuntimeError:
                self.shutdown_requested.set()
                break
            except Exception as queue_consumer_error:
                if self.shutdown_requested.is_set():
                    break
                time.sleep(0.1)

        self.status = 'Done'

    def create_plots(self, keys):
        self.plotdata = {key: np.array([], dtype=object) for key in keys}
        # figure out which key will be the x-axis
        if 'Step' not in keys:
            raise ValueError('The reporter has not step information, so there is no x-axis to plot graphs!')

    def remaining_time_reporter(self, data):
        self.classic_md_remain_data.append(data)

    def update_plot(self, msg):

        if self.plots_created and type(msg) == dict:
            for k, v in msg.items():
                current = self.plotdata.get(k)
                if current is None:
                    current = np.array([], dtype=object)

                updated = np.concatenate((current, np.array([v], dtype=object)), axis=None)
                if len(updated) > self.max_plot_points:
                    updated = updated[-self.max_plot_points:]

                if k == 'Progress (%)':
                    self.plotdata.update({k: updated})

                elif k == 'Reference MD Progress (%)':
                    try:
                        self.plotdata.pop('Progress (%)')
                    except:
                        pass
                    self.plotdata = {**{k: updated}, **self.plotdata}

                elif k == 'Dissipation MD Progress (%)':
                    try:
                        self.plotdata.pop('Reference MD Progress (%)')
                    except:
                        pass
                    self.plotdata = {**{k: updated}, **self.plotdata}
                else:
                    self.plotdata.update({k: updated})

            current_time = time.time()
            if current_time - self.last_plot_emit_time >= self.plot_emit_interval_sec:
                self._safe_emit(self.Signals.dataSignal, self.plotdata)
                self.last_plot_emit_time = current_time
            self._safe_emit(self.Signals.real_time_pertub_info, self.pymol_statu)

        if not self.plots_created and type(msg) == dict:
            self.create_plots(msg.keys())
            self.plots_created = True

        if type(msg) == list:
            self._safe_emit(self.Signals.decomp_process, msg)

        if type(msg) == str:
            if msg != "Progress Finished Succesfully :)":
                if not self._safe_emit(self.Signals.inform_about_situation, msg):
                    return
                if "error" in msg.lower() or "traceback" in msg.lower() or "exception" in msg.lower():
                    self._safe_emit(self.Signals.finish_alert, msg)
            else:
                self.decomp_data.append(int(100))
                # self.Signals.decomp_process.emit(self.decomp_data)
                self.update_plot(self.decomp_data)
                self._safe_emit(self.Signals.finish_alert, msg)
                self.plots_created = False


class Graphs(QWidget):
    TEMPERATURE_PEN = pg.mkPen((255, 0, 0), width=2)
    ENERGY_PEN = [pg.mkPen((255, 0, 0), width=2), pg.mkPen((0, 255, 0), width=2), pg.mkPen((0, 0, 255), width=2)]
    SPEED_PEN = pg.mkPen((200, 200, 200), width=2)
    TIME_PEN = pg.mkPen((255, 255, 0), width=2)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.real_time_as_minute = []
        self.real_speed = []
        self.current_step_keeper = None
        self.elapsed_time = 0.0

        self.setup_ui()

    def setup_ui(self):
        self.initialize_graph_layout()
        self.setup_temperature_graph()
        self.setup_energy_graph()
        # self.setup_simulation_speed_graph()
        # self.setup_simulation_time_graph()

    def initialize_graph_layout(self):
        self.configure_graph_options()
        self.win = pg.GraphicsLayoutWidget(show=False, title="RMSD Plotting Region")
        self.win.setStyleSheet("border: 2px solid green; padding: -5px; border-radius: 10px;")
        self.win.setWindowTitle('Real Time Simulation Monitoring')

    def configure_graph_options(self):
        pg.setConfigOption('background', None)
        pg.setConfigOption('foreground', (197, 198, 199))
        pg.setConfigOptions(antialias=True)

    def setup_temperature_graph(self, reset=None):
        if reset is None:
            self.temperature_graph = self.win.addPlot(title="Temperature")
            self.configure_graph(self.temperature_graph, 'Temperature', 'K')

        self.temperature_graph_plot = pg.PlotDataItem(clear=True, pen=self.TEMPERATURE_PEN,
                                                      name="Temperature", fillLevel=0.0, brush=(150, 150, 50, 30))
        self.temperature_graph.addItem(self.temperature_graph_plot)

    def setup_energy_graph(self, reset=None):
        if reset is None:
            self.energy_graph = self.win.addPlot(title="Energy")
            self.configure_graph(self.energy_graph, 'Energy', 'kJ/mole')

        self.potential_energy_graph = pg.PlotDataItem(clear=True, pen=self.ENERGY_PEN[0],
                                                      name="Potential", fillLevel=0.0, brush=(150, 150, 50, 30))
        self.kinetic_energy_graph = pg.PlotDataItem(clear=True, pen=self.ENERGY_PEN[1],
                                                    name="Kinetic", fillLevel=0.0, brush=(150, 150, 50, 30))
        self.total_energy_graph = pg.PlotDataItem(clear=True, pen=self.ENERGY_PEN[2],
                                                  name="Total", fillLevel=0.0, brush=(150, 150, 50, 30))
        self.energy_graph.addItem(self.potential_energy_graph)
        self.energy_graph.addItem(self.kinetic_energy_graph)
        self.energy_graph.addItem(self.total_energy_graph)
        # self.win.nextRow()

    def setup_simulation_speed_graph(self):
        self.simulation_speed_graph = self.win.addPlot(title="Speed", row=1, colspan=1)
        self.configure_graph(self.simulation_speed_graph, 'Speed', 'ns/day')
        self.simulation_speed_graph_plot = pg.PlotDataItem(clear=True, name="Speed (ns/day)", fillLevel=0.0,
                                                           brush=(150, 150, 50, 30), pen=self.SPEED_PEN,
                                                           symbolBrush=(255, 0, 0), symbolPen='w')
        self.simulation_speed_graph.addItem(self.simulation_speed_graph_plot)
        # self.win.nextRow()

    def setup_simulation_time_graph(self):
        self.simulation_time_graph = self.win.addPlot(title="Remaining Time", row=1, colspan=1)
        self.configure_graph(self.simulation_time_graph, 'Remaining Time', 'sec')
        self.simulation_time_graph_plot = pg.PlotDataItem(clear=True, pen=self.TIME_PEN,
                                                          fillLevel=0.0, name="Remaining Time (sec)",
                                                          brush=(150, 150, 50, 10))
        self.simulation_time_graph.addItem(self.simulation_time_graph_plot)
        # self.win.nextRow()

    def configure_graph(self, graph, title, units):
        graph.addLegend()
        graph.getViewBox().setBackgroundColor((129, 105, 161, 20))
        graph.setLabel('left', title, units=units)
        graph.setLabel('bottom', "Step")
        graph.setLogMode(x=False, y=False)
        graph.showGrid(x=True, y=True)
        graph.getAxis('left').enableAutoSIPrefix(False)

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

    def update_graph(self, data):

        try:
            x = np.array(data["Step"], dtype=np.float64)
            y_temp = np.array(data["Temperature (K)"], dtype=np.float64)
            y_potential = np.array(data["Potential Energy (kJ/mole)"], dtype=np.float64)
            y_kinetic = np.array(data["Kinetic Energy (kJ/mole)"], dtype=np.float64)
            y_total = np.array(data["Total Energy (kJ/mole)"], dtype=np.float64)

            # y_speed = np.array(data["Speed (ns/day)"])[-1]
            # real_speed = self.pretty_speed(y_speed)

            # time_remaining = np.array(data["Time Remaining"])[-1]
            # print(list(data.keys()))
            self.runner.Signals.run_speed.emit(data["Speed (ns/day)"])

            if list(data.keys())[0] == 'Progress (%)':
                self.runner.Signals.classic_md_remain_time.emit(data["Time Remaining"])

            if list(data.keys())[0] == 'Reference MD Progress (%)':
                self.runner.Signals.reference_md_remain_time.emit(data["Time Remaining"])

            if list(data.keys())[0] == 'Dissipation MD Progress (%)':
                self.runner.Signals.dissipation_md_remain_time.emit(data["Time Remaining"])

                if data["Time Remaining"][-1] == "0:00":
                    self.runner.Signals.run_speed.emit("--")

            if x.shape == y_temp.shape:
                if self.current_step_keeper is not None and len(self.current_step_keeper) > 1:
                    if self.current_step_keeper[-1] > x[-1]:
                        # if self.current_step_keeper[-1] - self.current_step_keeper[-2] > 1:
                        #     ticks = pg.VTickGroup(xvals=[self.current_step_keeper[-1]], yrange=[0, 2.5],
                        #                           pen={'color': 'g', 'width': 2.0, 'style': QtCore.Qt.DashLine})
                        #     self.temperature_graph_plot.getViewBox().addItem(ticks)
                        # x = np.append(self.current_step_keeper, self.current_step_keeper[-1] + 1)
                        length = max(len(x), len(self.current_step_keeper))
                        array1 = np.pad(x, (0, length - len(x)))
                        array2 = np.pad(self.current_step_keeper, (0, length - len(self.current_step_keeper)))
                        # Combine arrays using numpy.where
                        x = np.where(array2 > array1, array2, array1)
                        x[-1] = self.current_step_keeper[-1] + 1

                try:
                    # Ensure x and y arrays have the same shape before setting data
                    min_length = min(len(x), len(y_temp))
                    x = x[:min_length]
                    y_temp = y_temp[:min_length]
                    y_potential = y_potential[:min_length]
                    y_kinetic = y_kinetic[:min_length]
                    y_total = y_total[:min_length]

                    self.temperature_graph_plot.setData(x=x, y=y_temp)
                    self.temperature_graph.autoRange()
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
            else:
                pass

        except Exception as e:
            print("ERROR 2: %s" % e)

    def updating_decomposion(self, data_decomp):
        pass

    def updating_current_speed(self, data_speed):
        pass

    def run_script(self, contents):
        self.contents = contents
        self.runner = OpenMMScriptRunner(self.contents)
        self.runner.Signals.dataSignal.connect(lambda plotdata: self.update_graph(plotdata))
        # self.runner.Signals.decomp_process.connect(lambda decomp_data: self.updating_decomposion(decomp_data))
        # self.runner.Signals.run_speed.connect(lambda speed_data: self.updating_current_speed(speed_data))

    def stop_th(self):
        try:
            was_stopped = self.runner.stop_threads()
            if not was_stopped:
                return False

            self.runner.plotdata.clear()

            self.real_time_as_minute = []
            self.real_speed = []
            self.current_step_keeper = None
            self.elapsed_time = 0.0
            return True

        except Exception as Run_Stop_Error:
            Message_Boxes.Information_message(self, "There is no an active run!", str(Run_Stop_Error),
                                              Style.MessageBox_stylesheet)
            return False
# if __name__ == '__main__':
#     app = QtWidgets.QApplication(sys.argv)
#     main = Graphs(contents)
#     main.show()
#     sys.exit(app.exec_())
