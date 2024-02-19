from gui.ui_styles import Style
from src.mplwidget import *
from src.omm_runner import *
from src.pyside_dynamic import loadUi
import src.ui_functions as UIF

from PySide2 import QtCore
from PySide2.QtCore import Slot, QThread
from PySide2.QtWidgets import QMessageBox
import os
import pystache
from io import StringIO
import time
import queue
import threading
import itertools
import tokenize
from PySide2.QtCore import Signal
from ui_main import *


class Advanced(QtCore.QThread):

    def __init__(self, parent=None):
        super(Advanced, self).__init__(parent)
        self.start_monitoring = False

    @Slot()
    def send_arg_to_Engine(self):
        # self.stackedWidget.setCurrentIndex(1)
        global platform_name
        pdb_pfile = os.path.abspath(self.upload_pdb_lineEdit.text().strip()).replace('\\', '/')

        rigid_water = True
        minimize = True
        equilubrate = True
        no_minimize_value = True
        DCD_Reporter = True
        XTC_Reporter = False
        State_Data_Reporter = True
        Nonbounded_cutoff_active = True
        use_switching_distance = True
        Additional_Integrator = False
        properties_active = False
        CPU_properties_active = False
        Device_ID_active = False
        precision = 'single'
        water_active = False
        equilubrate_steps = self.Max_equilubrate_steps_textEdit.toPlainText()
        StateData_freq = int(self.StateData_frequency_lineEdit.text())
        # global selected_platform

        print("Simulation parameters preparing for the start ...")
        ## FORCEFIELD CONFIGURATIONS
        self.protein_ff = self.protein_forcefield_comboBox.currentText() + '.xml'

        if 'obc' or 'gbvi' in self.protein_forcefield_comboBox.currentText():
            water_active = True

        if self.protein_forcefield_comboBox.currentText() == 'charmm36':
            self.water_ff = '%s/%s.xml' % (
                self.protein_forcefield_comboBox.currentText(), self.water_forcefield_comboBox.currentText())
        else:
            self.water_ff = '%s.xml' % self.water_forcefield_comboBox.currentText()

        ## SOLUTION WATER MODEL
        if len((self.water_ff.split('/')[-1]).split('.')[0]) <= 5:
            self.water_model = (self.water_ff.split('/')[-1]).split('.')[0]

        elif (self.water_ff.split('/')[-1]).split('.')[0] == 'tip4pew':
            self.water_model = (self.water_ff.split('/')[-1]).split('.')[0]

        else:
            self.water_model = (self.water_ff.split('/')[-1]).split('.')[0][0:5]

        try:
            self.out_dir = self.Output_Folder_textEdit.toPlainText()
            if self.out_dir == '':
                """Output Directory Help."""
                msgbox = QMessageBox(QMessageBox.Question, "There is no specified output directory",
                                     "If you want to specify, please click the 'Yes' button, "
                                     "otherwise the system will create an output folder named 'output'")

                msgbox.setIcon(QMessageBox.Question)
                msgbox.addButton(QMessageBox.Yes)
                msgbox.addButton(QMessageBox.No)
                msgbox.setDefaultButton(QMessageBox.Yes)
                # msgbox.setCheckBox(cb)
                # msg.setWindowIcon(QIcon(ICON_PATH))
                msgbox.setStyleSheet(Style.MessageBox_stylesheet)
                msgbox.setWindowFlags(QtCore.Qt.FramelessWindowHint | QtCore.Qt.WindowStaysOnTopHint)
                output_answer = msgbox.exec_()

                if output_answer == QMessageBox.Yes:
                    return False

                if output_answer == QMessageBox.No:
                    from pathlib import Path
                    path_out = os.getcwd() + "/output"
                    Path(path_out).mkdir(parents=True, exist_ok=True)

                    self.out_dir = os.path.abspath(path_out.strip()).replace('\\', '/')
                    self.Output_Folder_textEdit.setText(self.out_dir)
            else:
                self.out_dir = os.path.abspath(self.out_dir.strip()).replace('\\', '/')
                self.Output_Folder_textEdit.setText(self.out_dir)
        except:
            QMessageBox.critical(self, "Error", "An error occured while getting output directory")

        if not self.rigid_water_checkBox.isChecked():
            rigid_water = False
        if not self.use_switching_checkBox.isChecked():
            use_switching_distance = False

        try:
            if use_switching_distance:
                if float(self.switching_distance_lineEdit.text().split('*')[0]) >= float(
                        self.nonbounded_CutOff_lineEdit.text().split('*')[0]):

                    """Display help."""
                    msg = QMessageBox()
                    # msg.setWindowTitle("Nonbonded Force Warning")
                    msg.setIcon(QMessageBox.Information)
                    msg.setText('Nonbonded Force Warning')
                    msg.setInformativeText("Switching Distance must smaller than nonbonded cutoff. "
                                           "Switching distance must satisfy \n0 <= r_switch < r_cutoff")

                    # msg.setWindowIcon(QIcon(ICON_PATH))
                    msg.setDetailedText("Solutions pages still in progress, so we can not redirect you a link.")
                    msg.setStyleSheet(Style.MessageBox_stylesheet)
                    msg.setWindowFlags(QtCore.Qt.FramelessWindowHint | Qt.WindowSystemMenuHint)
                    answer = msg.exec_()

                    if answer == QMessageBox.Ok:
                        return False

        except:
            QMessageBox.critical(self, "Error", "An error occured while getting Nonbonded Parameters")

        try:
            if self.Max_minimize_iter_lineEdit.text() == "":
                no_minimize_value = False

            if not self.minimize_checkBox.isChecked():
                minimize = False

            if not self.equilubrate_checkBox.isChecked():
                equilubrate = False

            if equilubrate_steps == "":
                equilubrate_steps = 500

            if not self.DCD_Reporter_checkBox.isChecked():
                DCD_Reporter = False

            if self.XTC_Reporter_checkBox.isChecked():
                XTC_Reporter = True

            if not self.State_Data_Reporter_checkBox.isChecked():
                State_Data_Reporter = False

            ## STATE DATA FREQUENCY
            if StateData_freq == "" and State_Data_Reporter:
                StateData_freq = 100
                if StateData_freq > self.Number_of_steps_spinBox.value():
                    StateData_freq = int(self.Number_of_steps_spinBox.value() / 2)
                    print("You did not enter the stateData frequency. Therefore it will be automatically set to %s."
                          % StateData_freq)

            elif StateData_freq != "" and State_Data_Reporter:
                if StateData_freq > self.Number_of_steps_spinBox.value():
                    StateData_freq = int(int(self.Number_of_steps_spinBox.value()) / 2)
                    print("You entered the stateData frequency bigger than total steps. Therefore it will be "
                          "automatically set to half of the total steps which is equal to %s. "
                          % StateData_freq)

            ## NONBOUNDED CUTOFF ACTIVE OR NOT
            if self.nonBounded_Method_comboBox.currentText == 'NoCutoff':
                Nonbounded_cutoff_active = False

            if self.Additional_Integrators_checkBox.isChecked():
                Additional_Integrator = True

            if self.equ_platform_comboBox.currentText() in ["CUDA", "OpenCL"]:
                properties_active = True

            if self.equ_platform_comboBox.currentText() == "CPU":
                CPU_properties_active = True

            if self.Device_ID_checkBox.isChecked():
                Device_ID_active = True

            platform_name, properties, precision = Advanced_Helper_Functions.selected_platform(self,
                                                                                               self.equ_platform_comboBox.currentText(),
                                                                                               Device_ID_active,
                                                                                               precision)
            cuda_precision_prefix = None
            cuda_active = False
            if platform_name == 'CUDA':
                cuda_precision_prefix = 'Cuda'
                cuda_active = True

        except Exception as e:
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
            print(exc_type, fname, exc_tb.tb_lineno)

        print("Parameters have been sent to OMM-Runner...")
        script_structure = dict(pdb=pdb_pfile,
                                output_folder=self.Output_Folder_textEdit.toPlainText(),
                                long_simulation_time=self.Number_of_steps_spinBox.value(),
                                long_simulation_time_unit=self.long_simulation_time_unit.currentText(),

                                Number_of_CPU=self.Number_CPU_spinBox_2.value(),
                                all_cpu=self.All_CPU_checkBox.isChecked(),
                                platform=platform_name, properties=properties, precision=precision,
                                cuda_active=cuda_active,
                                cuda_precision_prefix=cuda_precision_prefix,
                                properties_active=properties_active, CPU_properties_active=CPU_properties_active,
                                Device_ID_active=Device_ID_active,
                                Device_Number=self.Device_Number_comboBox.currentText(),

                                protein_ff=self.protein_ff,  # FORCEFIELD
                                water_ff=self.water_ff,  # FORCEFIELD
                                model_water=self.water_model,
                                water_active=water_active,  # FORCEFIELD if protein include obc or gbvi
                                water_padding=self.water_padding_lineEdit.text(),
                                perturbed_res_list=[self.selected_residues_listWidget.item(x).text()[:-1] for x in
                                                    range(self.selected_residues_listWidget.count())],
                                speed_factor=[int(r) if r.strip().isdigit() else r for r in
                                              self.R_factor_ComboBox.currentText()],
                                perturb_simulation_time=self.run_duration_spinBox.value(),
                                perturb_simulation_time_unit=self.perturb_time_unit_comboBox.currentText(),

                                integrator_kind=self.integrator_kind_comboBox.currentText(),  # INTEGRATOR
                                integrator_time_step=self.integrator_time_step_lineEdit.text().strip(),  # INTEGRATOR
                                integrator_time_step_unit=self.integrator_time_step_unit.currentText(),  # INTEGRATOR
                                friction=self.friction_lineEdit.text().strip(),  # ADDITIONAL INTEGRATOR
                                Temperature=self.temperature_lineEdit.text().strip(),  # ADDITIONAL INTEGRATOR
                                Additional_Integrator=Additional_Integrator,

                                NonBoundedMethod=self.nonBounded_Method_comboBox.currentText(),
                                Constraints=self.system_constraints_comboBox.currentText(),
                                Rigid_Water=rigid_water,
                                NonBounded_cutoff=self.nonbounded_CutOff_lineEdit.text(),
                                Nonbounded_cutoff_active=Nonbounded_cutoff_active,
                                use_switching_distance=use_switching_distance,
                                switching_distance=self.switching_distance_lineEdit.text(),
                                solvent_ionic_strength=self.ionic_strength_lineEdit.text(),

                                Number_of_steps=self.Number_of_steps_spinBox.value(),
                                Minimize=minimize,
                                no_minimize_value=no_minimize_value,
                                Max_minimization_iteration=self.Max_minimize_iter_lineEdit.text(),
                                Equilubrate=equilubrate, Equilubrate_steps=equilubrate_steps,
                                DCDReporter=DCD_Reporter, XTCReporter=XTC_Reporter,
                                DCD_write_freq=self.DCD_write_freq_lineEdit.text(),
                                DCD_output_name=self.DCD_Output_Name_lineEdit.text(),
                                XTC_write_freq=self.XTC_write_freq_lineEdit.text(),
                                XTC_output_name=self.XTC_Output_Name_lineEdit.text(),
                                State_Data_Reporter=State_Data_Reporter,
                                StateData_freq=StateData_freq,
                                output_directory=self.out_dir,
                                save_script=self.Save_Script_checkBox.isChecked(),
                                script_save_directory=self.script_name_lineEdit.text()
                                )

        self.created_script = Advanced_Helper_Functions.update_display(self, script_structure)

        self.start_monitoring = True
        return self.start_monitoring
        # return self.created_script
        # # self.update_display(script_structure)
        # self.Real_Time_Graphs.run_script(self.created_script)


class Advanced_Helper_Functions(QtCore.QThread):
    def __init__(self, parent=None):
        super(Advanced_Helper_Functions, self).__init__(parent)

    @Slot()
    def selected_platform(self, platform_name, Device_ID_active, precision):

        if platform_name == 'OpenCL' and Device_ID_active == True:
            properties = {'OpenCLPrecision': '%s' % precision,
                          'OpenCLDeviceIndex': '%s' % self.Device_Number_comboBox.currentText()}
            return platform_name, properties, precision

        if platform_name == 'OpenCL' and Device_ID_active == False:
            properties = {'OpenCLPrecision': '%s' % precision}
            return platform_name, properties, precision

        if platform_name == 'CUDA' and Device_ID_active == True:
            properties = {'CudaPrecision': '%s' % precision, 'CudaDeviceIndex': '%s' % self.Device_Index_Number}
            return platform_name, properties, precision

        if platform_name == 'CUDA' and not Device_ID_active:
            properties = {'CudaPrecision': '%s' % precision}
            return platform_name, properties, precision

        if platform_name == 'CPU':
            print("The CPU platform always uses 'mixed' precision.")
            print("Simulation process will use %s Thread(s)" % self.Number_CPU_spinBox_2.value())
            properties = {'CpuThreads': '%s' % self.Number_CPU_spinBox_2.value()}
            precision = 'mixed'
            return platform_name, properties, precision

        if platform_name == 'Reference':
            print("The Reference platform always uses 'double' precision.")
            properties = None
            precision = 'double'
            return platform_name, properties, precision
        print("System will use '%s' Platform with '%s' Precision" % (platform_name, precision))

    @Slot()
    def update_display(self, script_structure):
        renderer = pystache.Renderer()

        template = pystache.parse(open('src/template_script.txt').read())
        self.contents = renderer.render(template, script_structure)

        with open('readme.txt', 'w') as f:
            f.write(self.contents)

        if script_structure.get('save_script'):
            save_filename = script_structure.get('script_save_directory')
            # Check if the file already exists in the specified directory
            file_counter = 1
            while os.path.exists(os.path.join('temp', save_filename)):
                # If the file already exists, generate a new filename with a suffix
                save_filename, extension = os.path.splitext(save_filename)
                save_filename = f"{save_filename}_{file_counter}{extension}"
                file_counter += 1

            # Create or open a Python script file and write the content
            with open(os.path.join('temp', save_filename), 'w') as py_file:
                py_file.write(self.contents)

        return self.contents

