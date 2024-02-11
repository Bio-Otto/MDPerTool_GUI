from src.builder import *

from importlib import resources
import os
import sys
import platform
from PySide2.QtCore import (QCoreApplication, QPropertyAnimation, QDate, QDateTime, QMetaObject, QObject, QPoint, QRect,
                            QSize, QTime, QUrl, Qt, QEvent, QRegExp, QThreadPool, Signal)
from PySide2.QtGui import (QBrush, QColor, QConicalGradient, QCursor, QFont, QFontDatabase, QIcon, QKeySequence,
                           QLinearGradient, QPalette, QPainter, QPixmap, QRadialGradient, QIntValidator,
                           QRegExpValidator)
from PySide2 import QtWidgets, QtGui

# =================== > IMPORTS < =================== #
from matplotlib.backends.backend_qt5agg import (NavigationToolbar2QT as NavigationToolbar)
from openmm.app import Modeller
from pdbfixer import PDBFixer
import multiprocessing as mp

#  ==> GLOBALS
counter = 0


def center_window(widget):
    window = widget.window()
    window.setGeometry(
        QtWidgets.QStyle.alignedRect(
            QtCore.Qt.LeftToRight,
            QtCore.Qt.AlignCenter,
            window.size(),
            QtGui.QGuiApplication.primaryScreen().availableGeometry(),
        ),
    )


class PlotSignal(QObject):
    plot_network = Signal()


class MainWindow(QtWidgets.QMainWindow):
    global active_workers, network_holder, log_holder

    def __init__(self, parent=None):
        super(MainWindow, self).__init__(parent=parent)
        self.dragPos = None
        self.ui = loadUi('gui/MAIN_GUI.ui', self)

        ################################################################################################################
        #                                       ==> START OF WINDOW ATTRIBUTES <==                                     #
        ################################################################################################################
        # print('System: ' + platform.system())
        # print('Version: ' + platform.release())
        self.r_factor_count = 0
        # ------------------------------------ > START OF SIMULATION MONITORING < ------------------------------------ #
        self.created_script = None
        self.Real_Time_Graphs = Graphs()
        self.verticalLayout_22.addWidget(self.Real_Time_Graphs.win)
        # ------------------------------------- > END OF SIMULATION MONITORING < ------------------------------------- #

        # ############################################################################################################ #
        # -------------------------------------- > START OF MATPLOTLIB WIDGET < -------------------------------------- #
        # self.matplotlib_widget()
        self.plot_signal = PlotSignal()
        self.plot_signal.plot_network.connect(lambda: UIF.Functions.plot_networks(self))
        # -------------------------------------- > END OF MATPLOTLIB WIDGET < ---------------------------------------- #
        # ############################################################################################################ #

        # -------------------------------------- > Remove Standart Title Bar < --------------------------------------- #
        UIF.UIFunctions.removeTitleBar(True)
        # self.setWindowFlags(QtCore.Qt.FramelessWindowHint | QtCore.Qt.WindowStaysOnTopHint)

        # ----- > Set Window Title
        self.setWindowTitle('MDPerTool v0.1')
        # self.setWindowIcon(QIcon('%s/icons/MDPerTool.ico' % os.getcwd()))
        UIF.UIFunctions.labelTitle(self, 'MDPerTool - %s %s' % (platform.system(), platform.release()))
        UIF.UIFunctions.labelDescription(self, str(os.getcwd()))

        # # ----- > Move The Screen To Center
        # qtRectangle = self.frameGeometry()
        # centerPoint = QDesktopWidget().availableGeometry().center()
        # qtRectangle.moveCenter(centerPoint)
        # self.move(qtRectangle.topLeft())

        # ----- > Start With Standart Size
        startSize = QSize(1600, 1000)
        self.resize(startSize)
        self.setMinimumSize(startSize)
        # UIFunctions.enableMaximumSize(self, 500, 720)

        # --------------------------------------------- > CREATE MENUS < --------------------------------------------- #
        # ----- > Toggle Menu Size
        self.btn_toggle_menu.clicked.connect(lambda: UIF.UIFunctions.toggleMenu(self, 220, True))

        # ----- > Add Custom Menus
        self.stackedWidget.setMinimumWidth(20)
        UIF.UIFunctions.addNewMenu(self, "Perturbation", "btn_perturbation",
                                   "url(:/20x20/icons/20x20/chemical_20x20.png)", True)
        UIF.UIFunctions.addNewMenu(self, "Monitoring", "btn_monitoring",
                                   "url(:/20x20/icons/20x20/cil-monitor.png)", True)
        UIF.UIFunctions.addNewMenu(self, "Analysis", "btn_analysis",
                                   "url(:/16x16/icons/16x16/cil-chart-line.png)", True)
        """
        UIF.UIFunctions.addNewMenu(self, "Settings", "btn_settings",
                                   "url(:/20x20/icons/20x20/cil-equalizer.png)", False)
        """
        UIF.UIFunctions.addNewMenu(self, "About & Contact", "btn_about", "url(:/20x20/icons/20x20/cil-tag.png)", False)

        # ----- > Start Menu Selection
        UIF.UIFunctions.selectStandardMenu(self, "btn_perturbation")

        # ----- > Start Page
        self.stackedWidget.setCurrentWidget(self.page_home)

        # ----- > Show / Hide User Icon
        UIF.UIFunctions.userIcon(self, "HIO", "", True)

        # self.frame_label_top_btns.mouseMoveEvent = self.moveWindow

        # ----- > Load Definitions
        UIF.UIFunctions.uiDefinitions(self)

        ################################################################################################################
        #                                        ==> END OF WINDOW ATTRIBUTES <==                                      #
        ################################################################################################################

        # ------------------------------- > MAINWINDOW WIDGETS FUNCTIONS/PARAMETERS < -------------------------------- #
        self.quit_pushButton.clicked.connect(lambda: UIF.UIFunctions.close_application(self))
        self.PDB_ID_lineEdit.returnPressed.connect(self.fetch_and_load_pdbfile)
        self.stop_pushButton.clicked.connect(self.stop_button_clicked)
        self.upload_pdb_Button.clicked.connect(lambda: self.upload_pdb_from_local(manuel=True))
        self.Browse_Output_button.clicked.connect(lambda: self.output_folder_browse())
        self.PDB_ID_lineEdit.textChanged.connect(lambda: UIF.Functions.PDB_ID_lineEdit(self))
        self.fetch_pdb_Button.clicked.connect(lambda: self.fetch_and_load_pdbfile())
        self.integrator_kind_comboBox.currentTextChanged.connect(self.Stocasthic_Changed)
        UIF.Functions.Send_Available_Platforms_to_GUI(self)
        UIF.Functions.maximum_thread_of_system(self)
        self.mutator_checkBox.stateChanged.connect(lambda: UIF.Functions.mutator_isActive(self))
        self.platform_specific_precision_applying()

        self.node_threshold_checkBox.stateChanged.connect(lambda: UIF.Functions.node_threshold_use(self))
        self.equ_platform_comboBox.currentTextChanged.connect(
            lambda: UIF.Functions.platform_comboBox_Changed(self, eq_md=True, per_md=False))
        self.per_platform_comboBox.currentTextChanged.connect(
            lambda: UIF.Functions.platform_comboBox_Changed(self, eq_md=False, per_md=True))

        self.minimize_checkBox.stateChanged.connect(lambda: UIF.Functions.minimize_Step_isVisible(self))
        self.State_Data_Reporter_checkBox.stateChanged.connect(lambda: UIF.Functions.State_Data_Reporter_Changed(self))
        self.DCD_Reporter_checkBox.stateChanged.connect(lambda: UIF.Functions.DCD_Reporter_Changed(self))
        self.equilubrate_checkBox.stateChanged.connect(lambda: UIF.Functions.Equilibration_On_Off_Changed(self))
        self.Device_ID_checkBox.stateChanged.connect(lambda: UIF.Functions.UseDeviceID_On_Off_Changed(self))
        self.XTC_Reporter_checkBox.stateChanged.connect(lambda: UIF.Functions.XTC_Reporter_Changed(self))
        self.Run.clicked.connect(self.run_btn_clicked)
        self.load_sim_sample_pushButton.clicked.connect(lambda: UIF.Functions.load_sample_for_simulation(self))
        self.export_workspace_pushButton.clicked.connect(lambda: UIF.Functions.export_workspace(self))
        self.import_workspace_pushButton.clicked.connect(lambda: UIF.Functions.import_workspace(self))

        # --> RUN TIME SETTINGS
        self.run_duration_doubleSpinBox.valueChanged.connect(
            lambda: UIF.Functions.number_of_steps_changed_from_quick(self))
        self.long_simulation_time_unit.currentTextChanged.connect(
            lambda: UIF.Functions.number_of_steps_changed_from_quick(self))
        self.integrator_time_step.textChanged.connect(lambda: UIF.Functions.number_of_steps_changed_from_quick(self))
        self.Number_of_steps_spinBox.valueChanged.connect(
            lambda: UIF.Functions.number_of_steps_changed_from_advanced(self))

        # ============================================================================================================ #
        self.run = OpenMMScriptRunner

        self.run.Signals.decomp_process.connect(lambda decomp_data: self.progressBar_decomp.setValue(decomp_data[-1]))
        self.run.Signals.classic_md_remain_time.connect(
            lambda classic_md_remain_data: self.update_classic_md_remaining_time(classic_md_remain_data[-1]))
        self.run.Signals.reference_md_remain_time.connect(
            lambda reference_md_remain_data: self.update_reference_md_remaining_time(reference_md_remain_data[-1]))
        self.run.Signals.dissipation_md_remain_time.connect(
            lambda dissipation_md_remain_data: self.update_dissipation_md_remaining_time(dissipation_md_remain_data[-1]))

        self.run.Signals.run_speed.connect(lambda speed_data: self.update_simulation_speed(speed_data[-1]))


        self.run.Signals.finish_alert.connect(lambda finish_signal: self.finish_message(finish_signal))
        self.run.Signals.inform_about_situation.connect(
            lambda inform_message: self.inform_about_progress(inform_message))
        # ============================================================================================================ #

        # ------------------------------ > START OF ANALYSIS WINDOW RELEATED BUTTONS < ------------------------------- #
        self.response_time_upload_Button.clicked.connect(lambda: UIF.Functions.browse_responseTimeFile(self))
        self.response_time_lineEdit.textChanged.connect(self.response_time_graph_path_changed)
        self.source_res_comboBox.currentTextChanged.connect(self.response_time_graph_path_changed)
        self.all_targets_checkBox.stateChanged.connect(lambda: UIF.Functions.All_Residues_as_target_Changed(self))
        self.output_directory_button.clicked.connect(lambda: UIF.Functions.analysis_output_directory(self))
        self.upload_boundForm_pdb_Button.clicked.connect(lambda: self.upload_boundForm_pdb_from_local(manuel=True))
        self.add_residue_to_targets_pushButton.clicked.connect(lambda: UIF.Functions.add_residue_to_target_List(self))
        self.discard_residue_from_targets_pushButton.clicked.connect(
            lambda: UIF.Functions.discard_residue_from_target_List(self))
        self.get_conserv_score_pushButton.clicked.connect(lambda: UIF.Functions.get_conservation_scores(self))
        # self.show_2d_network_pushButton.clicked.connect(lambda: UIFunctions.start_VisJS_2D_Network(self))

        self.load_analysis_sample_pushButton.clicked.connect(lambda: UIF.Functions.load_sample_for_analysis(self))

        self.network_calculate_pushButton.clicked.connect(self.run_network_analysis)

        # ----------------------------------- > START OF PYMOL RELEATED BUTTONS < ------------------------------------ #
        UIF.UIFunctions.start_pymol(self)
        # UIFunctions.start_VisJS_2D_Network(self)
        self.add_residue_pushButton.clicked.connect(lambda: UIF.Functions.add_residue_toList(self))
        self.discard_residue_pushButton.clicked.connect(lambda: UIF.Functions.discard_residue_fromList(self))
        self.selected_residues_listWidget.itemDoubleClicked.connect(lambda: UIF.UIFunctions.show_residue_labels(self))
        self.refresh_pushButton.clicked.connect(lambda: UIF.UIFunctions.clear_residue_labels(self))
        self.activate_pymol_navigation.clicked.connect(lambda: UIF.UIFunctions.activate_navigation_on_Pymol(self))
        self.deactivate_pymol_navigation.clicked.connect(lambda: UIF.UIFunctions.deactivate_navigation_on_Pymol(self))
        self.visualization_Handel_buttons_changing()
        self.Handel_Buttons()
        self.ss_beatiful_snapshoot.clicked.connect(lambda: UIF.UIFunctions.show_beatiful_in_Pymol(self))
        self.get_figure_pushButton.clicked.connect(lambda: UIF.UIFunctions.save_as_png_Pymol(self))
        self.Handel_Save_Figure_Options_Changed()
        self.Handel_Save_Figure_Options()

        # ----> Analysis Menu Parameters
        self.threadpool = QThreadPool()
        self.active_workers = 0
        self.network_holder = []
        self.log_holder = []
        self.node_threshold = None
        self.thread_main = QtCore.QThread()
        self.thread_main.start()
        self.analysis_TabWidget.setTabsClosable(True)
        self.analysis_TabWidget.tabCloseRequested.connect(self.closeTab)

        # ################################################# TRAIL #################################################### #
        # self.add_tabb.clicked.connect(self.add_tab)
        self.analysis_TabWidget.tabBar().setTabButton(0, QTabBar.RightSide, None)
        self.analysis_TabWidget.tabBarClicked.connect(self.handle_tabbar_clicked)
        self.tab_count_on_analysis = str(self.analysis_TabWidget.count())

    def closeTab(self, currentIndex):
        self.analysis_TabWidget.removeTab(currentIndex)
        self.tab_count_on_analysis = self.analysis_TabWidget.count()

    def handle_tabbar_clicked(self, index_of_tab):
        pass

    # #################################################### TRAIL ##################################################### #
    def run_btn_clicked(self):
        self.r_factor_count = 0
        self.Run.setEnabled(False)
        self.start_monitoring = False
        self.start_monitoring = Advanced.send_arg_to_Engine(self)
        self.run.plotdata = {}

        if self.start_monitoring:
            self.show_simulation_monitoring()

        if not self.start_monitoring:
            self.Run.setEnabled(True)

    def platform_specific_precision_applying(self):
        eq_md_precission_indexes = {'single': self.eq_precision_comboBox.findText('single'),
                                    'mixed': self.eq_precision_comboBox.findText('mixed'),
                                    'double': self.eq_precision_comboBox.findText('double')}

        per_md_precission_indexes = {'single': self.per_precision_comboBox.findText('single'),
                                     'mixed': self.per_precision_comboBox.findText('mixed'),
                                     'double': self.per_precision_comboBox.findText('double')}

        self.equ_platform_comboBox.currentTextChanged.connect(
            lambda: UIF.Functions.precision_combobox_settings(self, eq_md_indexes=eq_md_precission_indexes,
                                                              per_md_indexes=None))
        self.per_platform_comboBox.currentTextChanged.connect(
            lambda: UIF.Functions.precision_combobox_settings(self, eq_md_indexes=None,
                                                              per_md_indexes=per_md_precission_indexes))
        UIF.Functions.precision_combobox_settings(self, eq_md_indexes=eq_md_precission_indexes,
                                                  per_md_indexes=per_md_precission_indexes)

    def stop_button_clicked(self):
        self.__stop = True
        try:
            self.Real_Time_Graphs.stop_th()
            self.Run.setEnabled(True)

            self.Real_Time_Graphs.temperature_graph_plot.clear()
            self.Real_Time_Graphs.simulation_speed_graph_plot.clear()

            self.Real_Time_Graphs.potential_energy_graph.clear()
            self.Real_Time_Graphs.kinetic_energy_graph.clear()
            self.Real_Time_Graphs.total_energy_graph.clear()

        except Exception as ins:
            pass
            # Message_Boxes.Information_message(self, "Info", str(ins), Style.MessageBox_stylesheet)

    def show_simulation_monitoring(self):
        # self.stackedWidget.setCurrentIndex(1)
        self.stackedWidget.setCurrentWidget(self.Simulation_Graphs)

        UIF.UIFunctions.resetStyle(self, "btn_monitoring")
        UIF.UIFunctions.labelPage(self, "MONITORING THE PROCESS")
        widget = self.findChild(QPushButton, "btn_monitoring")
        widget.setStyleSheet(UIF.UIFunctions.selectMenu(widget.styleSheet()))

        self.Real_Time_Graphs.run_script(self.created_script)

    def fetch_and_load_pdbfile(self):
        """
        Fetch action for getting crystal structure from pdb databank
        :return: path of downloaded pdb file or None conditions will return
        """
        fetched_and_modified_pdb = UIF.Functions.Fetch_PDB_File(self)

        if fetched_and_modified_pdb:
            UIF.UIFunctions.load_pdb_to_pymol(self, fetched_and_modified_pdb)

        if fetched_and_modified_pdb is None:
            UIF.Message_Boxes.Information_message(self, 'Wrong pdb id "%s"' % self.PDB_ID_lineEdit.text(),
                                                  'There is no protein crystal structure in this expression.',
                                                  Style.MessageBox_stylesheet)

    def upload_pdb_from_local(self, manuel):
        global selected_chains
        try:
            if manuel:
                upload_condition, pdb_path = UIF.Functions.browse_pdbFile(self)
            if not manuel:
                pdb_path = self.upload_pdb_lineEdit.text()

            if os.path.exists(pdb_path):
                fixer = PDBFixer(pdb_path)
                fixer.removeHeterogens(keepWater=False)

                modeller = Modeller(fixer.topology, fixer.positions)
                chains = [r.id for r in modeller.topology.chains()]

                checked_list = UIF.ChecklistDialog('Select the chain (s) to be used in the system', chains,
                                                   checked=True)
                pdb_fix_dialog_answer = checked_list.exec_()
                if pdb_fix_dialog_answer == QtWidgets.QDialog.Accepted:
                    selected_chains = [str(s) for s in checked_list.choices]

                    delete_chains = list(set(chains) - set(selected_chains))

                    modified_pdb = UIF.pdb_Tools.fetched_pdb_fix(self, pdb_path,
                                                                 self.Output_Folder_textEdit.toPlainText(),
                                                                 ph=7, chains_to_remove=delete_chains)

                    self.upload_pdb_lineEdit.setText(modified_pdb)

                    self.combobox = UIF.Helper_Functions.fill_residue_combobox(self, modified_pdb)
                    for i in self.combobox:
                        self.res1_comboBox.addItem(str(i))
                    self.res1_comboBox.clear()  # delete all items from comboBox
                    self.res1_comboBox.addItems(self.combobox)  # add the actual content of self.comboData

                elif pdb_fix_dialog_answer == QtWidgets.QDialog.Rejected:
                    modified_pdb = UIF.pdb_Tools.fetched_pdb_fix(self, pdb_path,
                                                                 self.Output_Folder_textEdit.toPlainText(),
                                                                 ph=7, chains_to_remove=None)

                    self.upload_pdb_textEdit.setText(modified_pdb)

                    self.combobox = UIF.Helper_Functions.fill_residue_combobox(self, modified_pdb)
                    for i in self.combobox:
                        self.res1_comboBox.addItem(str(i))
                    self.res1_comboBox.clear()  # delete all items from comboBox
                    self.res1_comboBox.addItems(self.combobox)  # add the actual content of self.comboData

                UIF.UIFunctions.load_pdb_to_pymol(self, modified_pdb)

        except Exception as Error:
            print("ERROR: ", Error)
            pass

    def upload_boundForm_pdb_from_local(self, manuel=False):
        global selected_chains
        try:
            if manuel:
                upload_condition, pdb_path = UIF.Functions.browse_bound_form_pdbFile(self)
            if not manuel:
                pdb_path = self.boundForm_pdb_lineedit.text()

            if os.path.exists(pdb_path):
                fixer = PDBFixer(pdb_path)
                fixer.removeHeterogens(keepWater=False)

                modeller = Modeller(fixer.topology, fixer.positions)
                chains = [r.id for r in modeller.topology.chains()]

                checked_list = UIF.ChecklistDialog('Select the chain (s) to be used in the system', chains,
                                                   checked=True)
                pdb_fix_dialog_answer = checked_list.exec_()
                if pdb_fix_dialog_answer == QtWidgets.QDialog.Accepted:
                    selected_chains = [str(s) for s in checked_list.choices]

                    delete_chains = list(set(chains) - set(selected_chains))

                    modified_pdb = UIF.pdb_Tools.fetched_pdb_fix(self, pdb_path,
                                                                 self.Output_Folder_textEdit.toPlainText(),
                                                                 ph=7, chains_to_remove=delete_chains)

                    self.boundForm_pdb_lineedit.setText(modified_pdb)

                    self.target_combobox = UIF.Helper_Functions.fill_residue_combobox(self, modified_pdb)

                    for i in self.target_combobox:
                        self.target_res_comboBox.addItem(str(i))
                        self.source_res_comboBox.addItem(str(i))

                    self.target_res_comboBox.clear()  # delete all items from comboBox
                    self.target_res_comboBox.addItems(self.target_combobox)  # add the actual content of self.comboData
                    self.source_res_comboBox.clear()  # delete all items from comboBox
                    self.source_res_comboBox.addItems(self.target_combobox)  # add the actual content of self.comboData

                elif pdb_fix_dialog_answer == QtWidgets.QDialog.Rejected:
                    modified_pdb = UIF.pdb_Tools.fetched_pdb_fix(self, pdb_path,
                                                                 self.Output_Folder_textEdit.toPlainText(),
                                                                 ph=7, chains_to_remove=None)

                    self.boundForm_pdb_lineedit.setText(modified_pdb)

                    self.target_combobox = UIF.Helper_Functions.fill_residue_combobox(self, modified_pdb)
                    for i in self.target_combobox:
                        self.target_res_comboBox.addItem(str(i))
                        self.source_res_comboBox.addItem(str(i))

                    self.target_res_comboBox.clear()  # delete all items from comboBox
                    self.target_res_comboBox.addItems(self.target_combobox)  # add the actual content of self.comboData
                    self.source_res_comboBox.clear()  # delete all items from comboBox
                    self.source_res_comboBox.addItems(self.target_combobox)  # add the actual content of self.comboData

                UIF.UIFunctions.load_pdb_to_3DNetwork(self, modified_pdb)

        except TypeError:
            pass

    def output_folder_browse(self):
        UIF.Functions.output_file(self)

    def Stocasthic_Changed(self):
        UIF.Functions.Stochastic_changed(self)

    def finish_message(self, alert_message):
        self.r_factor_count += 1
        if self.r_factor_count == len(str(self.R_factor_ComboBox.currentText())):
            UIF.Message_Boxes.Succesfully_message(self, "Thumbs Up :)", alert_message, Style.MessageBox_stylesheet)
            self.Run.setEnabled(True)

    def inform_about_progress(self, message):
        message_regular = '<font color="yellow">%s</font>' % message
        self.label_top_info_1.setText(message_regular)

    def update_classic_md_remaining_time(self, time):
        font = QFont("Segoe UI", 10, QFont.Bold)  # Segoe UI, 10 boyut, bold font
        font.setLetterSpacing(QFont.AbsoluteSpacing, 1)  # Dikey hizalamayı ayarla
        remainingTime_regular = '<font color="red"><b>%s</b></font>' % time
        self.label_56.setText(remainingTime_regular)
        self.label_56.setFont(font)
        self.label_56.setAlignment(Qt.AlignVCenter)  # Dikey hizalamayı ayarla

    def update_reference_md_remaining_time(self, time):
        font = QFont("Segoe UI", 10, QFont.Bold)  # Segoe UI, 10 boyut, bold font
        font.setLetterSpacing(QFont.AbsoluteSpacing, 1)  # Dikey hizalamayı ayarla
        remainingTime_regular = '<font color="red"><b>%s</b></font>' % time
        self.label_57.setText(remainingTime_regular)
        self.label_57.setFont(font)
        self.label_57.setAlignment(Qt.AlignVCenter)  # Dikey hizalamayı ayarla

    def update_dissipation_md_remaining_time(self, time):
        font = QFont("Segoe UI", 10, QFont.Bold)  # Segoe UI, 10 boyut, bold font
        font.setLetterSpacing(QFont.AbsoluteSpacing, 1)  # Dikey hizalamayı ayarla
        remainingTime_regular = '<font color="red"><b>%s</b></font>' % time
        self.label_58.setText(remainingTime_regular)
        self.label_58.setFont(font)
        self.label_58.setAlignment(Qt.AlignVCenter)  # Dikey hizalamayı ayarla

    def update_simulation_speed(self, speed):
        try:
            speed_float = float(speed)
            formatted_speed = '{:.2f} ns/day'.format(speed_float)
        except ValueError:
            formatted_speed = '-- ns/day'

        font = QFont("Segoe UI", 10, QFont.Bold)
        font.setLetterSpacing(QFont.AbsoluteSpacing, 1)
        run_realtime_speed = '<font color="red"><b>{}</b></font>'.format(formatted_speed)

        self.label_59.setText(run_realtime_speed)
        self.label_59.setFont(font)
        self.label_59.setAlignment(Qt.AlignVCenter)

    ####################################################################################################################
    #                                     ==> START OF DYNAMIC MENUS FUNCTIONS < ==                                    #
    ####################################################################################################################
    def Button(self):
        # GET BT CLICKED
        btnWidget = self.sender()

        # PAGE HOME
        if btnWidget.objectName() == "btn_home":
            self.stackedWidget.setCurrentWidget(self.page_home)
            UIF.UIFunctions.resetStyle(self, "btn_home")
            UIF.UIFunctions.labelPage(self, "Home")
            btnWidget.setStyleSheet(UIF.UIFunctions.selectMenu(btnWidget.styleSheet()))

        # PAGE PERTURBATION
        if btnWidget.objectName() == "btn_perturbation":
            self.stackedWidget.setCurrentWidget(self.page_home)
            UIF.UIFunctions.resetStyle(self, "btn_perturbation")
            UIF.UIFunctions.labelPage(self, "Perturbation User Interface")
            btnWidget.setStyleSheet(UIF.UIFunctions.selectMenu(btnWidget.styleSheet()))

        # PAGE MONITORING
        if btnWidget.objectName() == "btn_monitoring":
            self.stackedWidget.setCurrentWidget(self.Simulation_Graphs)
            UIF.UIFunctions.resetStyle(self, "btn_monitoring")
            UIF.UIFunctions.labelPage(self, "MONITORING THE PROCESS")
            btnWidget.setStyleSheet(UIF.UIFunctions.selectMenu(btnWidget.styleSheet()))

        # PAGE ANALYSIS
        if btnWidget.objectName() == "btn_analysis":
            self.stackedWidget.setCurrentWidget(self.Analysis)
            UIF.UIFunctions.resetStyle(self, "btn_analysis")
            UIF.UIFunctions.labelPage(self, "ENERGY DISSIPATION ANALYSIS")
            btnWidget.setStyleSheet(UIF.UIFunctions.selectMenu(btnWidget.styleSheet()))

        # PAGE WIDGETS
        if btnWidget.objectName() == "btn_settings":
            self.stackedWidget.setCurrentWidget(self.page_settings)
            UIF.UIFunctions.resetStyle(self, "btn_settings")
            UIF.UIFunctions.labelPage(self, "Settings")
            btnWidget.setStyleSheet(UIF.UIFunctions.selectMenu(btnWidget.styleSheet()))

        # PAGE ABOUT
        if btnWidget.objectName() == "btn_about":
            self.stackedWidget.setCurrentWidget(self.page_contact)
            UIF.UIFunctions.resetStyle(self, "btn_about")
            UIF.UIFunctions.labelPage(self, "ABOUT & CONTACT")
            btnWidget.setStyleSheet(UIF.UIFunctions.selectMenu(btnWidget.styleSheet()))

    ####################################################################################################################
    #                                     == > END OF DYNAMIC MENUS FUNCTIONS < ==                                     #
    ####################################################################################################################
    ####################################################################################################################
    #                                          == > START OF APP EVENTS < ==                                           #
    ####################################################################################################################

    # def moveEvent(self, e):
    #     super(MainWindow, self).moveEvent(e)

    # ----- > Move Window / Maximize / Restore
    @staticmethod
    def moveWindow(self, event):
        # ----- > If Maximized Change to Normal
        if UIF.UIFunctions.returStatus() == 1:
            UIF.UIFunctions.maximize_restore(self)
        # ----- > Move Window
        if event.buttons() == Qt.LeftButton:
            self.move(self.pos() + event.globalPos() - self.dragPos)
            self.dragPos = event.globalPos()
            event.accept()

    # ----- > Mouse Click Event
    def mousePressEvent(self, event):
        self.dragPos = event.globalPos()
        self.old_width = self.width()
        self.old_height = self.height()
        """
        if event.buttons() == Qt.LeftButton:
            print('Mouse click: LEFT CLICK')
        if event.buttons() == Qt.RightButton:
            print('Mouse click: RIGHT CLICK')
        if event.buttons() == Qt.MidButton:
            print('Mouse click: MIDDLE BUTTON')
        """

    # ---> Mouse Click Event - Start
    def mouseMoveEvent(self, event):
        if self.dragPos:
            delta = QPoint(event.globalPos() - self.dragPos)
            if (self.dragPos.x() > self.x() + self.old_width - 10) or (
                    self.dragPos.y() > self.y() + self.old_height - 10):
                self.setFixedSize(self.old_width + delta.x(), self.old_height + delta.y())
            else:
                self.move(self.x() + delta.x(), self.y() + delta.y())
                self.dragPos = event.globalPos()

        if UIF.UIFunctions.returStatus() == 1:
            UIF.UIFunctions.maximize_restore(self)

    def mouseReleaseEvent(self, event):
        self.dragPos = None

    # ----- > Mouse Double Click
    def mouseDoubleClickEvent(self, event):
        widget = self.childAt(event.pos())
        if widget.objectName() == 'label_title_bar_top':
            if self.isMaximized():
                self.showNormal()
                self.btn_maximize_restore.setToolTip("Maximize")
                self.btn_maximize_restore.setIcon(QtGui.QIcon(u":/16x16/icons/16x16/cil-window-maximize.png"))
                self.frame_size_grip.show()
            else:
                self.showMaximized()
                self.horizontalLayout.setContentsMargins(0, 0, 0, 0)
                self.btn_maximize_restore.setToolTip("Restore")
                self.btn_maximize_restore.setIcon(QtGui.QIcon(u":/16x16/icons/16x16/cil-window-restore.png"))
                self.frame_size_grip.hide()

    # ----- > Key Pressed
    def keyPressEvent(self, event):
        # if event.key() == Qt.Key_Escape:
        #     self.close()
        # print('Key: ' + str(event.key()) + ' | Text Press: ' + str(event.text()))
        pass

    # ----- > Resize Event
    def resizeEvent(self, event):
        self.resizeFunction()
        path = QtGui.QPainterPath()
        path.addRoundedRect(QtCore.QRectF(self.rect()), 5, 5)
        reg = QtGui.QRegion(path.toFillPolygon().toPolygon())
        self.setMask(reg)
        return super(MainWindow, self).resizeEvent(event)

    # ----- > When Resize The Frame It's Height and Width Will Print
    def resizeFunction(self):
        # print('Height: ' + str(self.height()) + ' | Width: ' + str(self.width()))
        pass

    ##################################################################################################################
    #                                           == > END OF APP EVENTS < ==                                          #
    ##################################################################################################################
    ##################################################################################################################
    #                                   == > START OF PYMOL NAVIGATION TOOLBAR < ==                                  #
    ##################################################################################################################
    def visualization_Handel_buttons_changing(self):
        self.hide_visualization_settings()

    def Handel_Buttons(self):
        self.show_navigation.clicked.connect(self.show_visualization_settings)
        self.hide_navigation.clicked.connect(self.hide_visualization_settings)

    def show_visualization_settings(self):
        self.visualization_settings_groupBox.show()
        self.show_navigation.hide()
        self.hide_navigation.show()

    def hide_visualization_settings(self):
        self.visualization_settings_groupBox.hide()
        self.show_navigation.show()
        self.hide_navigation.hide()

    # ----------------------------------------- > FIGURE OPTIONS IN PYMOL < ------------------------------------------ #
    def Handel_Save_Figure_Options_Changed(self):
        self.hide_figure_options()

    def Handel_Save_Figure_Options(self):
        self.save_as_png_pushButton.clicked.connect(self.show_figure_options)
        self.hide_figure_settings_pushButton.clicked.connect(self.hide_figure_options)
        self.width_horizontalSlider.valueChanged.connect(self.figure_width_label)
        self.height_horizontalSlider.valueChanged.connect(self.figure_height_label)
        self.dpi_horizontalSlider.valueChanged.connect(self.figure_dpi_label)
        self.ray_horizontalSlider.valueChanged.connect(self.figure_ray_label)

    def show_figure_options(self):
        self.figure_settings_groupBox.show()

    def hide_figure_options(self):
        self.figure_settings_groupBox.hide()

    def figure_width_label(self):
        self.pymol_width_label.setText("Width: " + str(self.width_horizontalSlider.value()))

    def figure_height_label(self):
        self.pymol_height_label.setText("Height: " + str(self.height_horizontalSlider.value()))

    def figure_dpi_label(self):
        self.pymol_dpi_label.setText("Dpi: " + str(self.dpi_horizontalSlider.value()))

    def figure_ray_label(self):
        self.pymol_ray_label.setText("Ray: " + str(self.ray_horizontalSlider.value()))

    ####################################################################################################################
    #                                   == > END OF PYMOL NAVIGATION TOOLBAR < ==                                      #
    ####################################################################################################################
    def matplotlib_widget(self):
        self.matplotlib_widget = WidgetPlot(self)
        self.verticalLayout_23.addWidget(self.matplotlib_widget.toolbar)
        self.verticalLayout_23.addWidget(self.matplotlib_widget.canvas)

    def response_time_graph_path_changed(self):
        possible_path = str(self.response_time_lineEdit.text())
        if os.path.exists(possible_path.strip()) and possible_path.split('.')[-1] == 'csv':
            source_residue = self.source_res_comboBox.currentText()
            row, col, Response_Count = getResponseTimeGraph(possible_path)
            """
            if source_residue == '':
                self.matplotlib_widget.canvas.plot(Response_Count, source_residue=None)
            if source_residue != '':
                self.matplotlib_widget.canvas.plot(Response_Count, source_residue=source_residue)
            """

    def run_network_analysis(self):

        # try:
        if os.path.exists(self.response_time_lineEdit.text()) and os.path.exists(self.boundForm_pdb_lineedit.text()):
            output_directory_for_network = self.net_output_directory_lineedit.text()
            if not os.path.exists(output_directory_for_network):
                """Output Directory Help."""
                msgbox = QMessageBox(QMessageBox.Question, "There is no specified output directory",
                                     "There is no specified output directory\n\n"
                                     "If you want to specify, please click the 'Yes' button, "
                                     "otherwise the system will create an output folder named 'out_for_net_analysis'")

                msgbox.setIcon(QMessageBox.Question)
                msgbox.addButton(QMessageBox.Yes)
                msgbox.addButton(QMessageBox.No)
                msgbox.setDefaultButton(QMessageBox.Yes)
                msgbox.setStyleSheet(Style.MessageBox_stylesheet)
                msgbox.setWindowFlags(QtCore.Qt.FramelessWindowHint | QtCore.Qt.WindowStaysOnTopHint)
                output_answer = msgbox.exec_()

                if output_answer == QMessageBox.Yes:
                    return False

                if output_answer == QMessageBox.No:
                    from pathlib import Path
                    path_out = os.getcwd() + "/out_for_net_analysis"
                    Path(path_out).mkdir(parents=True, exist_ok=True)

                    output_directory_for_network = os.path.abspath(path_out.strip()).replace('\\', '/')
                    self.net_output_directory_lineedit.setText(output_directory_for_network)
            else:
                output_directory_for_network = os.path.abspath(output_directory_for_network.strip()).replace('\\', '/')
                self.net_output_directory_lineedit.setText(output_directory_for_network)

            UIF.Functions.calculate_intersection_network(self)
        # except:
        #     exc_type, exc_obj, exc_tb = sys.exc_info()
        #     fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        #     print(exc_type, fname, exc_tb.tb_lineno)
        #     QMessageBox(QMessageBox.Critical, "Error", "Problem\nnAn error occured while getting output directory")

    def network_calc_thread(self, x):
        # print('mycallback is called with {}'.format(x))
        # progress = QProgressDialog("Downloading...", "Abort", 0, 0, parent=self)
        # progress.setWindowModality(QtCore.Qt.WindowModal)
        # progress.show()
        # time.sleep(0.1)
        # while thread_1.is_alive():
        #     QApplication([]).processEvents()
        #     # progress.setLabelText(thread.status())
        #     if progress.wasCanceled():
        #         sys.exit()

        # progress.setValue(total)
        # progress.setLabelText(thread.status())
        # QtGui.QMessageBox.information(self, "Done", "Download is complete.")
        # QMessageBox.information(self, "Done", "thread.status()")
        # progress.close()
        # return True
        pass

    def load_nx_to_VisJS_2D_Network(self, intersection_graph_html='2d_network.html', intersection_gml_file=None):
        initial_2d_network_html_directory = os.path.join(os.getcwd(), 'analysis')
        initial_2d_network_html_path = os.path.join(initial_2d_network_html_directory, intersection_graph_html)
        self.VisJSEngineView.load_network_component(network=intersection_gml_file,
                                                    html_file=initial_2d_network_html_path)
        self.VisJSEngineView()

    ####################################################################################################################
    #                                       == > START OF SPLASH SCREEN < ==                                           #
    ####################################################################################################################


class SplashScreen(QMainWindow):

    def __init__(self):
        QMainWindow.__init__(self)
        self.ui = loadUi('gui/splash_screen.ui', self)

        # ----- > Set App Icon
        app_icon = QtGui.QIcon()
        app_icon.addFile('icon.png')

        # ----- > Remove Background and Give icon to The Program
        self.setWindowIcon(QIcon(app_icon))
        self.setWindowFlags(QtCore.Qt.FramelessWindowHint | QtCore.Qt.WindowStaysOnTopHint)
        self.setAttribute(QtCore.Qt.WA_TranslucentBackground)

        # ----- > Drop Shadow Effect
        self.shadow = QGraphicsDropShadowEffect(self)
        self.shadow.setBlurRadius(20)
        self.shadow.setXOffset(0)
        self.shadow.setYOffset(0)
        self.shadow.setColor(QColor(0, 0, 0, 60))
        self.dropShadowFrame.setGraphicsEffect(self.shadow)

        # ----- > Qtimer For Change Description on Splash Screen
        self.timer = QtCore.QTimer()
        self.timer.timeout.connect(self.progress)
        self.timer.start(19)

        # ----- > Change Description on Splash Screen
        self.label_title.setText("MDPerTool v0.0.1")
        self.label_credits.setText("<strong>Devoloped</strong> by Ozbek's Lab")

        # ----- > Initial Text on Splash Screen
        self.label_description.setText("<strong>WELCOME</strong> TO MDPerTool V0.0.1 PLATFORM")

        # ----- > CHANGE TEXT ON SPLASH SCREEN
        QtCore.QTimer.singleShot(1500, lambda: self.label_description.setText("<strong>LOADING</strong> ENVIRONMENT"))
        QtCore.QTimer.singleShot(2000,
                                 lambda: self.label_description.setText("<strong>LOADING</strong> USER INTERFACE"))

        self.show()

    # ------------------------------------------ > START OF APP FUNCTIONS < ------------------------------------------ #
    def progress(self):
        """
            Splash Screen Progress Bar
        :return: There is no returm
        """
        global counter

        self.progressBar.setValue(counter)
        if counter > 100:
            self.timer.stop()
            self.main = MainWindow()
            QtCore.QTimer.singleShot(0, lambda: center_window(self.main))
            self.main.show()
            self.close()

        counter += 1

    # ------------------------------------------- > END OF APP FUNCTIONS < ------------------------------------------- #


def run_gui():
    mp.freeze_support()
    app = QtWidgets.QApplication(sys.argv)
    app.setStyleSheet(Style.QToolTip_stylesheet)
    window = SplashScreen()
    QtCore.QTimer.singleShot(0, lambda: center_window(window))
    sys.exit(app.exec_())


if __name__ == "__main__":
    if os.name == 'nt':
        # print("PLATFORM IS WINDOWS ..")
        import PySide2

        pyqt = os.path.dirname(PySide2.__file__)
        QApplication.addLibraryPath(os.path.join(pyqt, "plugins"))

    mp.freeze_support()
    app = QtWidgets.QApplication(sys.argv)
    app.setStyleSheet(Style.QToolTip_stylesheet)
    window = SplashScreen()
    QtCore.QTimer.singleShot(0, lambda: center_window(window))
    sys.exit(app.exec_())
