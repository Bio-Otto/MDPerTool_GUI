import os.path
import tempfile
from functools import partial
import networkx as nx
from PySide2.QtWidgets import (
    QFileDialog,
    QDialog,
    QTableWidgetItem,
    QTabBar,
    QVBoxLayout,
    QHBoxLayout,
    QListWidgetItem,
    QCheckBox,
    QGroupBox,
    QLabel,
    QSpinBox,
)
from PySide2.QtCore import Qt, Slot, QMutexLocker, QMutex
from PySide2.QtGui import QColor
from PySide2 import QtCore, QtGui, QtWidgets
from pathlib import Path
from pdbfixer import PDBFixer
from os import path
from urllib.request import urlretrieve
from openmm.app import PDBFile, Modeller
from gui.ui_styles import Style
from .checkBox_menu import ChecklistDialog
from .message import Message_Boxes
from .PyMolWidget import PymolQtWidget
from ._advanced_platform_options import (
    initialize_advanced_platform_options,
    get_advanced_platform_properties,
)
import multiprocessing as mp
from analysis.pdbsum_conservation_puller import get_conservation_scores
from analysis.createRNetwork import (Pymol_Visualize_Path, Shortest_Path_Visualize, get_residues)
from analysis.pathway_analysis import (
    summarize_target_pathways,
    build_critical_residue_rows,
    count_reachability,
    extract_target_graph_pairs,
    build_done_message,
)
from analysis.network_summary_service import build_network_summary_rows, build_union_graph, build_intersection_graph
from analysis.network_motif_service import build_motif_summary_rows
from analysis.network_significance_service import build_network_significance_rows
from analysis.response_dynamics_service import build_response_dynamics_payload
from analysis.network_workflow_service import (
    prepare_general_network_engine_from_ui,
    prepare_network_ui_for_run,
    start_general_network_background_worker,
    handle_general_network_ready_callback,
    present_network_failure_warning,
)
from analysis.analysis_presenters import (
    present_warning_payload,
    populate_pathway_summary_table,
    populate_critical_residue_table,
    populate_reachability_qc,
    populate_residue_response_table,
    populate_domain_summary_table,
    populate_network_summary_table,
    populate_motif_summary_table,
    populate_significance_table,
    populate_metrics_table,
    populate_qc_table,
    populate_provenance_table,
)
from .config import write_output_configuration_file, read_output_configuration_file, config_template
from .mplwidget import getResponseTimeGraph, WidgetPlot
from src.file_dialog import Dialog as file_dialog

# Import modularized helpers
from .helpers import (
    UILayoutManager,
    PyMOLVisualizer,
    ResidueManager,
    NetworkParametersManager,
    ProgressDialogManager,
)

# List of three-letter amino acid residue codes
amino_acid_residues = [
    "ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE",
    "LYS", "LEU", "MET", "ASN", "PRO", "GLN", "ARG", "SER",
    "THR", "VAL", "TRP", "TYR"
]

colors = ['#957DAD', '#D291BC', '#8565c4', '#8dbdc7', '#B3ABCF', '#b5b1c8', '#e8abb5', '#77dd77', '#aa9499',
          '#b39eb5', '#c6a4a4', '#ff694f', '#95b8d1', '#52b2cf', '#d3ab9e', '#fb6f92', '#872187',
          '#74138C', '', '#ff70a6', '#dab894', '#f6bc66', '#e27396', '#6e78ff', '#ff686b']


def _safe_getcwd():
    try:
        return os.getcwd()
    except FileNotFoundError:
        return str(Path.home())


class GeneralNetworkWorker(QtCore.QObject):
    finished = QtCore.Signal(object, object, int, object)
    failed = QtCore.Signal(str)

    def __init__(self, engine):
        super().__init__()
        self.engine = engine

    @QtCore.Slot()
    def run(self):
        try:
            initial_network, res_id_list, len_of_retimes = self.engine.calculate_general_network()
            self.finished.emit(initial_network, res_id_list, len_of_retimes, self.engine)
        except Exception as error:
            self.failed.emit(str(error))


class Helper_Functions():
    """
    Facade class that delegates to modularized helper classes.
    Maintains backward compatibility while providing cleaner separation of concerns.
    Supports both static and instance method usage patterns.
    """
    # Class-level helpers for static method calls
    _ui_manager_class = UILayoutManager()
    _pymol_viz_class = None  # Initialized on first use
    _residue_manager_class = ResidueManager()
    _network_manager_class = NetworkParametersManager()
    _progress_dialog_manager_class = None  # Initialized on first use with stylesheet

    def __init__(self, main_window=None):
        """Initialize helper instances for different concerns."""
        self._main_window = main_window
        self._ui_manager = UILayoutManager()
        self._pymol_viz = PyMOLVisualizer() if main_window else None
        self._residue_manager = ResidueManager()
        self._network_manager = NetworkParametersManager()
        self._progress_dialog_manager = None

    @staticmethod
    def _get_progress_dialog_manager(parent=None):
        """Get or create the class-level progress dialog manager."""
        if Helper_Functions._progress_dialog_manager_class is None:
            Helper_Functions._progress_dialog_manager_class = ProgressDialogManager(
                parent, Style.QProgressDialog_stylesheet
            )
        return Helper_Functions._progress_dialog_manager_class

    def fill_residue_combobox(self, pdb_path):
        from prody.proteins.pdbfile import parsePDB
        combobox = []

        pdb = parsePDB(pdb_path)
        protein = pdb.select('protein')
        for model in protein.getHierView():
            for chain in model:
                combobox.append(str(chain).replace(" ", "") + str(model).split(" ")[1])

        return combobox

    def available_platforms(self):
        """Delegate to NetworkParametersManager for platform detection."""
        return Helper_Functions._network_manager_class.available_platforms()

    def show_visualization_settings_on_analysis(self, analysis_settings_groupBox, show_navigation_button,
                                                hide_navigation_button):
        """Delegate to UILayoutManager."""
        Helper_Functions._ui_manager_class.show_visualization_settings_on_analysis(analysis_settings_groupBox, show_navigation_button,
                                                                                   hide_navigation_button)

    def hide_visualization_settings_on_analysis(self, analysis_settings_groupBox, show_navigation_button,
                                                hide_navigation_button):
        """Delegate to UILayoutManager."""
        Helper_Functions._ui_manager_class.hide_visualization_settings_on_analysis(show_navigation_button, hide_navigation_button,
                                                                                   analysis_settings_groupBox)

    def visualization_Handle_buttons_changing_on_analysis(self, analysis_settings_groupBox, show_navigation_button,
                                                          hide_navigation_button):
        """Delegate to UILayoutManager."""
        Helper_Functions._ui_manager_class.hide_visualization_settings_on_analysis(show_navigation_button, hide_navigation_button,
                                                                                   analysis_settings_groupBox)

    def Handle_Buttons_on_analysis(self, analysis_settings_groupBox, show_navigation_button, hide_navigation_button):
        """Delegate to UILayoutManager."""
        Helper_Functions._ui_manager_class.Handle_Buttons_on_analysis(analysis_settings_groupBox, show_navigation_button,
                                                                      hide_navigation_button)

    def activate_navigation_on_Pymol(self, created_pymol_widget):
        """Delegate to PyMOLVisualizer."""
        if not Helper_Functions._pymol_viz_class:
            Helper_Functions._pymol_viz_class = PyMOLVisualizer()
        Helper_Functions._pymol_viz_class.activate_navigation_on_Pymol(created_pymol_widget)

    def deactivate_navigation_on_Pymol(self, created_pymol_widget):
        """Delegate to PyMOLVisualizer."""
        if not Helper_Functions._pymol_viz_class:
            Helper_Functions._pymol_viz_class = PyMOLVisualizer()
        Helper_Functions._pymol_viz_class.deactivate_navigation_on_Pymol(created_pymol_widget)

    def clear_residue_labels(self, created_pymol_widget):
        """Delegate to PyMOLVisualizer."""
        if not Helper_Functions._pymol_viz_class:
            Helper_Functions._pymol_viz_class = PyMOLVisualizer()
        Helper_Functions._pymol_viz_class.clear_residue_labels(created_pymol_widget)

    def show_beautiful_in_Pymol(self, created_pymol_widget):
        """Delegate to PyMOLVisualizer."""
        if not Helper_Functions._pymol_viz_class:
            Helper_Functions._pymol_viz_class = PyMOLVisualizer()
        Helper_Functions._pymol_viz_class.show_beautiful_in_Pymol(created_pymol_widget)

    def save_as_png_Pymol(self, created_pymol_widget, width_horizontalSlider, height_horizontalSlider,
                          dpi_horizontalSlider, ray_horizontalSlider):
        """Delegate to PyMOLVisualizer (with main_window context for dialogs)."""
        if not Helper_Functions._pymol_viz_class:
            Helper_Functions._pymol_viz_class = PyMOLVisualizer()
        Helper_Functions._pymol_viz_class.save_as_png_Pymol(self, created_pymol_widget, width_horizontalSlider,
                                                            height_horizontalSlider, dpi_horizontalSlider, ray_horizontalSlider)

    # ----------------------------------------- > FIGURE OPTIONS IN PYMOL < ------------------------------------------ #
    def Handle_Save_Figure_Options_on_analysis_Changed(self, figure_settings_groupBox_on_analysis):
        """Delegate to UILayoutManager."""
        Helper_Functions._ui_manager_class.hide_figure_options_on_analysis(figure_settings_groupBox_on_analysis)

    def Handle_Save_Figure_Options_on_analysis(self, save_as_png_pushButton, hide_figure_settings_pushButton,
                                               width_horizontalSlider, height_horizontalSlider, dpi_horizontalSlider,
                                               ray_horizontalSlider, figure_settings_groupBox_on_analysis,
                                               pymol_width_label, pymol_height_label, pymol_dpi_label, pymol_ray_label):
        """Delegate to UILayoutManager."""
        Helper_Functions._ui_manager_class.Handle_Save_Figure_Options_on_analysis(
            save_as_png_pushButton, hide_figure_settings_pushButton,
            width_horizontalSlider, height_horizontalSlider, dpi_horizontalSlider,
            ray_horizontalSlider, figure_settings_groupBox_on_analysis,
            pymol_width_label, pymol_height_label, pymol_dpi_label, pymol_ray_label)

    def show_figure_options_on_analysis(self, figure_settings_groupBox_on_analysis):
        """Delegate to UILayoutManager."""
        Helper_Functions._ui_manager_class.show_figure_options_on_analysis(figure_settings_groupBox_on_analysis)

    def hide_figure_options_on_analysis(self, figure_settings_groupBox_on_analysis):
        """Delegate to UILayoutManager."""
        Helper_Functions._ui_manager_class.hide_figure_options_on_analysis(figure_settings_groupBox_on_analysis)

    def figure_width_label_on_analysis(self, width_horizontalSlider, with_label):
        """Delegate to UILayoutManager."""
        Helper_Functions._ui_manager_class.figure_width_label_on_analysis(width_horizontalSlider, with_label)

    def figure_height_label_on_analysis(self, height_horizontalSlider, height_label):
        """Delegate to UILayoutManager."""
        Helper_Functions._ui_manager_class.figure_height_label_on_analysis(height_horizontalSlider, height_label)

    def figure_dpi_label_on_analysis(self, dpi_horizontalSlider, dpi_label):
        """Delegate to UILayoutManager."""
        Helper_Functions._ui_manager_class.figure_dpi_label_on_analysis(dpi_horizontalSlider, dpi_label)

    def figure_ray_label_on_analysis(self, ray_horizontalSlider, ray_label):
        """Delegate to UILayoutManager."""
        Helper_Functions._ui_manager_class.figure_ray_label_on_analysis(ray_horizontalSlider, ray_label)

    def create_shortest_path_string(self, sp):
        shrotest_str_form = ''
        for res_id in range(len(sp)):
            if res_id == len(sp) - 1:
                shrotest_str_form += '%s' % sp[res_id]
            else:
                shrotest_str_form += '%s --> ' % sp[res_id]
        return shrotest_str_form

class Functions:
    global active_workers

    # ########################################### ANALYSIS WINDOW FUNCTIONS ############################################
    @Slot()
    def on_started(self):
        """Called when each worker thread starts."""
        return

    @Slot()
    def progress_fn(self, progress_on_net_calc):
        # print(progress_on_net_calc, "IS DONE")
        if self.active_workers == 0:
            pass

    def thread_complete(self):
        self._completed_workers = getattr(self, '_completed_workers', 0) + 1
        expected_workers = max(0, getattr(self, '_expected_workers', 0))

        if expected_workers > 0 and self._completed_workers >= expected_workers:
            if getattr(self, '_network_finalize_done', False):
                return
            self._network_finalize_done = True
            # Close progress dialog when all threads complete
            Helper_Functions._get_progress_dialog_manager().close(process_events=True)
            # Re-enable main window after calculation is done
            self.setEnabled(True)
            if hasattr(self, '_active_network_engine'):
                self._active_network_engine = None
            if len(self.network_holder) != 0:
                self.plot_signal.plot_network.emit()

            # Show PC computation results if they exist in the output directory
            try:
                out_dir = getattr(self, '_network_output_directory', getattr(self, 'output_directory', ''))
                if out_dir and os.path.exists(out_dir):
                    import glob
                    pc_plots = glob.glob(os.path.join(out_dir, "PC_Score_*_propagation_plot.png"))
                    if pc_plots:
                        from PySide2.QtWidgets import QMessageBox
                        msg = QMessageBox()
                        msg.setIcon(QMessageBox.Information)
                        msg.setWindowTitle("Network Metric Calculation")
                        msg.setStyleSheet(Style.MessageBox_stylesheet)
                        msg.setAttribute(QtCore.Qt.WA_StyledBackground, True)
                        msg.setText(f"Propagation Coefficient (PC) analysis effectively completed!\n\nFound {len(pc_plots)} plot(s) in:\n{out_dir}")
                        msg.addButton(QMessageBox.Ok)
                        btn_open = msg.addButton("Show Plots", QMessageBox.ActionRole)
                        msg.exec_()
                        
                        if msg.clickedButton() == btn_open:
                            import subprocess
                            for plot in pc_plots:
                                if os.name == 'nt':
                                    os.startfile(plot)
                                else:
                                    subprocess.call(['xdg-open', plot])
            except Exception as e:
                print(f"[Warning] Error opening PC Plots: {e}")

    def print_output(self, s):
        self.network_holder.append(s[0])
        self.log_holder.append(s[1])

    def calculate_intersection_network(self):
        """Start network calculation with Qt-safe background worker."""
        # Prepare progress dialog and UI state for background execution
        progress_manager = Helper_Functions._get_progress_dialog_manager(self)
        self.progress = prepare_network_ui_for_run(
            target=self,
            progress_manager=progress_manager,
        )

        try:
            engine = prepare_general_network_engine_from_ui(self, atom_type='CA')

            start_general_network_background_worker(
                target=self,
                engine=engine,
                thread_class=QtCore.QThread,
                worker_class=GeneralNetworkWorker,
                on_ready=partial(Functions._on_general_network_ready, self),
                on_failed=partial(Functions._on_general_network_failed, self),
            )

        except Exception as error:
            present_network_failure_warning(
                target=self,
                progress_manager=progress_manager,
                error_text=str(error),
                phase='startup',
                show_warning_fn=Message_Boxes.Warning_message,
                stylesheet=Style.MessageBox_stylesheet,
                warning_presenter=present_warning_payload,
            )

    def _on_general_network_ready(self, initial_network, resId_List, len_of_reTimes, engine):
        progress_manager = Helper_Functions._get_progress_dialog_manager()
        self.initial_network = initial_network

        warning_payload = handle_general_network_ready_callback(
            target=self,
            progress_manager=progress_manager,
            initial_network=self.initial_network,
            residue_ids=resId_List,
            response_time_count=len_of_reTimes,
            engine=engine,
            threadpool=self.threadpool,
            on_progress=partial(Functions.progress_fn, self),
            on_started=partial(Functions.on_started, self),
            on_result=partial(Functions.print_output, self),
            on_finished=partial(Functions.thread_complete, self),
            on_progress_message=self.inform_about_progress,
            output_directory=getattr(self, '_network_output_directory', self.output_directory),
        )
        if warning_payload is not None:
            present_warning_payload(
                show_warning_fn=Message_Boxes.Warning_message,
                owner=self,
                payload=warning_payload,
                stylesheet=Style.MessageBox_stylesheet,
            )
            return

    def _on_general_network_failed(self, error_text):
        progress_manager = Helper_Functions._get_progress_dialog_manager()
        present_network_failure_warning(
            target=self,
            progress_manager=progress_manager,
            error_text=str(error_text),
            phase='runtime',
            show_warning_fn=Message_Boxes.Warning_message,
            stylesheet=Style.MessageBox_stylesheet,
            warning_presenter=present_warning_payload,
        )

    def plot_networks(self):

        clean_graph_list = []
        clean_log_list = []

        if self.node_threshold is not None:
            # print("node threshold is NOT --> NONE")
            for k, i in enumerate(self.network_holder):
                if len(i.nodes()) > self.node_threshold:
                    clean_log_list.append(self.log_holder[k])
                    clean_graph_list.append(i)

        if self.node_threshold is None:
            # print("node threshold is --> NONE")

            for k, i in enumerate(self.network_holder):
                if len(i.nodes()) > 0:
                    clean_log_list.append(self.log_holder[k])
                    clean_graph_list.append(i)

        # CREATE AN INTERSECTION GRAPH AND WRITE TO GML FILE
        if len(clean_graph_list) > 0:
            all_graph_list = clean_graph_list
            tab = QtWidgets.QWidget()
            tab.setObjectName("Analysis_" + str(self.tab_count_on_analysis))

            self.analysis_TabWidget.tabBar().setTabButton(0, QTabBar.RightSide, None)
            self.tab_count_on_analysis = self.analysis_TabWidget.count()

            horizontalLayout = QtWidgets.QHBoxLayout(tab)
            horizontalLayout.setObjectName("horizontalLayout_" + str(self.tab_count_on_analysis))
            gridLayout = QtWidgets.QGridLayout()
            gridLayout.setObjectName("gridLayout_" + str(self.tab_count_on_analysis))

            analysis_data_tabwidget = QtWidgets.QTabWidget(tab)
            analysis_data_tabwidget.setObjectName("analysis_data_tabwidget")
            analysis_data_tabwidget.setMinimumSize(QtCore.QSize(450, 520))
            analysis_data_tabwidget.setMaximumSize(QtCore.QSize(520, 16777215))

            paths_tab = QtWidgets.QWidget()
            paths_tab.setObjectName("analysis_paths_tab")
            paths_tab_layout = QtWidgets.QVBoxLayout(paths_tab)
            paths_tab_layout.setContentsMargins(0, 0, 0, 0)
            paths_tab_layout.setObjectName("analysis_paths_tab_layout")

            label_2 = QtWidgets.QLabel(paths_tab)
            label_2.setText("Shortest Path(s)")
            label_2.setMinimumSize(QtCore.QSize(0, 22))
            label_2.setMaximumSize(QtCore.QSize(16777215, 22))
            label_2.setStyleSheet("QLabel {\n"
                                  "    background-color: rgb(27, 29, 35);\n"
                                  "    border-radius: 5px;\n"
                                  "    border: 2px solid rgb(27, 29, 35);\n"
                                  "    padding: 1px 1px 1px 1px;\n"
                                  "    \n"
                                  "    border-bottom-color: rgb(157, 90, 198);\n"
                                  "}\n"
                                  "\n"
                                  "\n"
                                  "QLabel:hover{\n"
                                  "    border: 2px solid rgb(64, 71, 88);\n"
                                  "    selection-color: rgb(127, 5, 64);\n"
                                  "\n"
                                  "}")
            label_2.setObjectName("label_29")
            paths_tab_layout.addWidget(label_2)

            shortest_path_listWidget = QtWidgets.QListWidget(paths_tab)
            shortest_path_listWidget.setMaximumSize(QtCore.QSize(16777215, 16777215))
            shortest_path_listWidget.setObjectName("shortest_path_listWidget")
            paths_tab_layout.addWidget(shortest_path_listWidget)

            tables_tab = QtWidgets.QWidget()
            tables_tab.setObjectName("analysis_tables_tab")
            tables_tab_layout = QtWidgets.QVBoxLayout(tables_tab)
            tables_tab_layout.setContentsMargins(0, 0, 0, 0)
            tables_tab_layout.setObjectName("analysis_tables_tab_layout")

            tables_tabwidget = QtWidgets.QTabWidget(tables_tab)
            tables_tabwidget.setObjectName("analysis_tables_tabwidget")
            tables_tab_layout.addWidget(tables_tabwidget)

            plot_tab = QtWidgets.QWidget()
            plot_tab.setObjectName("analysis_plot_tab")
            plot_tab_layout = QtWidgets.QVBoxLayout(plot_tab)
            plot_tab_layout.setContentsMargins(0, 0, 0, 0)
            plot_tab_layout.setObjectName("analysis_plot_tab_layout")

            metrics_tab = QtWidgets.QWidget()
            metrics_tab.setObjectName("analysis_metrics_tab")
            metrics_tab_layout = QtWidgets.QVBoxLayout(metrics_tab)
            metrics_tab_layout.setContentsMargins(0, 0, 0, 0)
            metrics_tab_layout.setObjectName("analysis_metrics_tab_layout")

            qc_tab = QtWidgets.QWidget()
            qc_tab.setObjectName("analysis_qc_tab")
            qc_tab_layout = QtWidgets.QVBoxLayout(qc_tab)
            qc_tab_layout.setContentsMargins(0, 0, 0, 0)
            qc_tab_layout.setObjectName("analysis_qc_tab_layout")

            provenance_tab = QtWidgets.QWidget()
            provenance_tab.setObjectName("analysis_provenance_tab")
            provenance_tab_layout = QtWidgets.QVBoxLayout(provenance_tab)
            provenance_tab_layout.setContentsMargins(0, 0, 0, 0)
            provenance_tab_layout.setObjectName("analysis_provenance_tab_layout")

            metrics_table = QtWidgets.QTableWidget(metrics_tab)
            metrics_table.setObjectName("response_metrics_table")
            metrics_table.setColumnCount(2)
            metrics_table.setRowCount(11)
            metrics_table.setHorizontalHeaderLabels(["Metric", "Value"])
            metrics_table.horizontalHeader().setStretchLastSection(True)
            metrics_table.verticalHeader().setVisible(False)
            metrics_table.setMaximumSize(QtCore.QSize(16777215, 16777215))
            metrics_table.setMinimumSize(QtCore.QSize(0, 0))
            metrics_table.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
            metrics_tab_layout.addWidget(metrics_table)

            qc_table = QtWidgets.QTableWidget(qc_tab)
            qc_table.setObjectName("response_qc_table")
            qc_table.setColumnCount(2)
            qc_table.setRowCount(6)
            qc_table.setHorizontalHeaderLabels(["QC Metric", "Value"])
            qc_table.horizontalHeader().setStretchLastSection(True)
            qc_table.verticalHeader().setVisible(False)
            qc_table.setMaximumSize(QtCore.QSize(16777215, 16777215))
            qc_table.setMinimumSize(QtCore.QSize(0, 0))
            qc_table.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
            qc_tab_layout.addWidget(qc_table)

            provenance_table = QtWidgets.QTableWidget(provenance_tab)
            provenance_table.setObjectName("analysis_provenance_table")
            provenance_table.setColumnCount(2)
            provenance_table.setRowCount(4)
            provenance_table.setHorizontalHeaderLabels(["Field", "Value"])
            provenance_table.horizontalHeader().setStretchLastSection(True)
            provenance_table.verticalHeader().setVisible(False)
            provenance_table.setMaximumSize(QtCore.QSize(16777215, 16777215))
            provenance_table.setMinimumSize(QtCore.QSize(0, 0))
            provenance_table.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
            provenance_tab_layout.addWidget(provenance_table)

            # --> Response Dynamics Tab
            response_dynamics_tab = QtWidgets.QWidget()
            response_dynamics_tab.setObjectName("analysis_response_dynamics_tab")
            response_dynamics_tab_layout = QtWidgets.QVBoxLayout(response_dynamics_tab)
            response_dynamics_tab_layout.setContentsMargins(0, 0, 0, 0)
            response_dynamics_tab_layout.setObjectName("analysis_response_dynamics_tab_layout")

            response_recalc_controls_layout = QtWidgets.QHBoxLayout()
            response_recalc_controls_layout.setObjectName("analysis_response_recalc_controls_layout")

            response_threshold_label = QtWidgets.QLabel(response_dynamics_tab)
            response_threshold_label.setObjectName("analysis_response_threshold_label")
            response_threshold_label.setText("Response Threshold:")

            response_threshold_spinbox = QtWidgets.QDoubleSpinBox(response_dynamics_tab)
            response_threshold_spinbox.setObjectName("analysis_response_threshold_spinbox")
            response_threshold_spinbox.setDecimals(3)
            response_threshold_spinbox.setRange(0.0, 10.0)
            response_threshold_spinbox.setSingleStep(0.005)
            response_threshold_spinbox.setValue(float(getattr(self, '_analysis_response_threshold', 0.01)))
            response_threshold_spinbox.setMinimumHeight(33)
            response_threshold_spinbox.setStyleSheet("""QDoubleSpinBox
            {
                color: white;
                border: 2px solid rgb(52, 59, 72);
                border-radius: 5px;
                background-color: rgb(27, 29, 35);
                border-width: 1px;
                padding: 4px 8px;
                selection-background-color: rgb(110, 105, 225);
            }
            QDoubleSpinBox:hover
            {
                border: 2px solid rgb(61, 70, 86);
            }
            QDoubleSpinBox:focus
            {
                border: 2px solid rgb(22, 200, 244);
            }
            QDoubleSpinBox::up-button, QDoubleSpinBox::down-button
            {
                width: 16px;
                border: none;
                background-color: rgb(52, 59, 72);
            }
            QDoubleSpinBox::up-button:hover, QDoubleSpinBox::down-button:hover
            {
                background-color: rgb(110, 105, 225);
            }""")

            recalculate_response_button = QtWidgets.QPushButton(response_dynamics_tab)
            recalculate_response_button.setObjectName("analysis_recalculate_response_button")
            recalculate_response_button.setText("Recalculate Response Time")
            recalculate_response_button.setMinimumSize(QtCore.QSize(190, 33))
            recalculate_response_button.setMaximumSize(QtCore.QSize(230, 33))
            recalculate_response_button.setStyleSheet("""QPushButton
            {
                color: white;
                border: 2px solid rgb(52, 59, 72);
                border-radius: 5px;
                background-color: rgb(253, 1, 136);
                border-width: 1px;
                padding: 5px;
                outline: none;
            }
            QPushButton:hover
            {
                background-color: rgb(22, 200, 244);
                border: 2px solid rgb(61, 70, 86);
            }
            QPushButton:pressed
            {
                background-color: rgb(15, 133, 163);
                border: 2px solid rgb(43, 50, 61);
            }""")
            recalculate_response_button.setToolTip("Recompute responseTimes from energy CSV files using selected threshold")

            response_recalc_controls_layout.addWidget(response_threshold_label)
            response_recalc_controls_layout.addWidget(response_threshold_spinbox)
            response_recalc_controls_layout.addWidget(recalculate_response_button)
            response_recalc_controls_layout.addStretch(1)
            response_dynamics_tab_layout.addLayout(response_recalc_controls_layout)

            response_dynamics_tabwidget = QtWidgets.QTabWidget(response_dynamics_tab)
            response_dynamics_tabwidget.setObjectName("analysis_response_dynamics_tabwidget")
            response_dynamics_tab_layout.addWidget(response_dynamics_tabwidget)

            # Per-residue response time table
            residue_response_tab = QtWidgets.QWidget()
            residue_response_tab.setObjectName("analysis_residue_response_tab")
            residue_response_tab_layout = QtWidgets.QVBoxLayout(residue_response_tab)
            residue_response_tab_layout.setContentsMargins(0, 0, 0, 0)
            residue_response_tab_layout.setObjectName("analysis_residue_response_tab_layout")

            residue_response_table = QtWidgets.QTableWidget(residue_response_tab)
            residue_response_table.setObjectName("residue_response_table")
            residue_response_table.setColumnCount(3)
            residue_response_table.setHorizontalHeaderLabels(["Residue ID", "Name", "Response Frame"])
            residue_response_table.horizontalHeader().setStretchLastSection(True)
            residue_response_table.verticalHeader().setVisible(False)
            residue_response_table.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
            residue_response_table.setMaximumSize(QtCore.QSize(16777215, 16777215))
            residue_response_table.setMinimumSize(QtCore.QSize(0, 280))
            residue_response_tab_layout.addWidget(residue_response_table)
            response_dynamics_tabwidget.addTab(residue_response_tab, "Per-Residue")

            # Domain summary table
            domain_summary_tab = QtWidgets.QWidget()
            domain_summary_tab.setObjectName("analysis_domain_summary_tab")
            domain_summary_tab_layout = QtWidgets.QVBoxLayout(domain_summary_tab)
            domain_summary_tab_layout.setContentsMargins(0, 0, 0, 0)
            domain_summary_tab_layout.setObjectName("analysis_domain_summary_tab_layout")

            domain_summary_table = QtWidgets.QTableWidget(domain_summary_tab)
            domain_summary_table.setObjectName("domain_summary_table")
            domain_summary_table.setColumnCount(5)
            domain_summary_table.setHorizontalHeaderLabels(["Domain", "# Residues", "Mean (ps)", "Min (ps)", "Max (ps)"])
            domain_summary_table.horizontalHeader().setStretchLastSection(True)
            domain_summary_table.verticalHeader().setVisible(False)
            domain_summary_table.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
            domain_summary_table.setMaximumSize(QtCore.QSize(16777215, 16777215))
            domain_summary_table.setMinimumSize(QtCore.QSize(0, 280))
            domain_summary_tab_layout.addWidget(domain_summary_table)
            response_dynamics_tabwidget.addTab(domain_summary_tab, "Domain Summary")

            # Target-level pathway summary table
            pathway_summary_tab = QtWidgets.QWidget()
            pathway_summary_tab.setObjectName("analysis_pathway_summary_tab")
            pathway_summary_tab_layout = QtWidgets.QVBoxLayout(pathway_summary_tab)
            pathway_summary_tab_layout.setContentsMargins(0, 0, 0, 0)
            pathway_summary_tab_layout.setObjectName("analysis_pathway_summary_tab_layout")

            pathway_summary_table = QtWidgets.QTableWidget(pathway_summary_tab)
            pathway_summary_table.setObjectName("pathway_summary_table")
            pathway_summary_table.setColumnCount(4)
            pathway_summary_table.setHorizontalHeaderLabels(["Target", "Nodes", "Shortest Path", "Route Type"])
            pathway_summary_table.horizontalHeader().setStretchLastSection(True)
            pathway_summary_table.verticalHeader().setVisible(False)
            pathway_summary_table.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
            pathway_summary_table.setMaximumSize(QtCore.QSize(16777215, 16777215))
            pathway_summary_table.setMinimumSize(QtCore.QSize(0, 280))
            pathway_summary_tab_layout.addWidget(pathway_summary_table)
            response_dynamics_tabwidget.addTab(pathway_summary_tab, "Pathway Analysis")

            # Critical residue ranking table
            critical_residue_tab = QtWidgets.QWidget()
            critical_residue_tab.setObjectName("analysis_critical_residue_tab")
            critical_residue_tab_layout = QtWidgets.QVBoxLayout(critical_residue_tab)
            critical_residue_tab_layout.setContentsMargins(0, 0, 0, 0)
            critical_residue_tab_layout.setObjectName("analysis_critical_residue_tab_layout")

            critical_residue_table = QtWidgets.QTableWidget(critical_residue_tab)
            critical_residue_table.setObjectName("critical_residue_table")
            critical_residue_table.setColumnCount(4)
            critical_residue_table.setHorizontalHeaderLabels(["Residue", "Path Hits", "Betweenness", "Composite Score"])
            critical_residue_table.horizontalHeader().setStretchLastSection(True)
            critical_residue_table.verticalHeader().setVisible(False)
            critical_residue_table.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
            critical_residue_table.setMaximumSize(QtCore.QSize(16777215, 16777215))
            critical_residue_table.setMinimumSize(QtCore.QSize(0, 280))
            critical_residue_tab_layout.addWidget(critical_residue_table)
            response_dynamics_tabwidget.addTab(critical_residue_tab, "Critical Residues")

            network_summary_tab = QtWidgets.QWidget()
            network_summary_tab.setObjectName("analysis_network_summary_tab")
            network_summary_tab_layout = QtWidgets.QVBoxLayout(network_summary_tab)
            network_summary_tab_layout.setContentsMargins(0, 0, 0, 0)
            network_summary_tab_layout.setObjectName("analysis_network_summary_tab_layout")

            network_summary_table = QtWidgets.QTableWidget(network_summary_tab)
            network_summary_table.setObjectName("network_summary_table")
            network_summary_table.setColumnCount(9)
            network_summary_table.setHorizontalHeaderLabels([
                "Network",
                "Nodes",
                "Edges",
                "Radius",
                "Diameter",
                "Characteristic Path",
                "Shortest Paths",
                "Avg Neighbors",
                "Clustering",
            ])
            network_summary_table.horizontalHeader().setStretchLastSection(True)
            network_summary_table.verticalHeader().setVisible(False)
            network_summary_table.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
            network_summary_table.setMaximumSize(QtCore.QSize(16777215, 16777215))
            network_summary_table.setMinimumSize(QtCore.QSize(0, 220))
            network_summary_tab_layout.addWidget(network_summary_table)
            response_dynamics_tabwidget.addTab(network_summary_tab, "Network Summary")

            motif_summary_tab = QtWidgets.QWidget()
            motif_summary_tab.setObjectName("analysis_motif_summary_tab")
            motif_summary_tab_layout = QtWidgets.QVBoxLayout(motif_summary_tab)
            motif_summary_tab_layout.setContentsMargins(0, 0, 0, 0)
            motif_summary_tab_layout.setObjectName("analysis_motif_summary_tab_layout")

            motif_summary_table = QtWidgets.QTableWidget(motif_summary_tab)
            motif_summary_table.setObjectName("motif_summary_table")
            motif_summary_table.setColumnCount(6)
            motif_summary_table.setHorizontalHeaderLabels([
                "Size",
                "Motif ID",
                "Edges",
                "Occurrence",
                "Frequency",
                "Scope",
            ])
            motif_summary_table.horizontalHeader().setStretchLastSection(True)
            motif_summary_table.verticalHeader().setVisible(False)
            motif_summary_table.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
            motif_summary_table.setMaximumSize(QtCore.QSize(16777215, 16777215))
            motif_summary_table.setMinimumSize(QtCore.QSize(0, 220))
            motif_summary_tab_layout.addWidget(motif_summary_table)
            response_dynamics_tabwidget.addTab(motif_summary_tab, "Motif Summary")

            significance_tab = QtWidgets.QWidget()
            significance_tab.setObjectName("analysis_significance_tab")
            significance_tab_layout = QtWidgets.QVBoxLayout(significance_tab)
            significance_tab_layout.setContentsMargins(0, 0, 0, 0)
            significance_tab_layout.setObjectName("analysis_significance_tab_layout")

            significance_table = QtWidgets.QTableWidget(significance_tab)
            significance_table.setObjectName("significance_table")
            significance_table.setColumnCount(8)
            significance_table.setHorizontalHeaderLabels([
                "Label",
                "Population",
                "Candidate Set",
                "Overlap",
                "P-value",
                "Candidate Nodes",
                "Overlap Nodes",
                "Description",
            ])
            significance_table.horizontalHeader().setStretchLastSection(True)
            significance_table.verticalHeader().setVisible(False)
            significance_table.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
            significance_table.setMaximumSize(QtCore.QSize(16777215, 16777215))
            significance_table.setMinimumSize(QtCore.QSize(0, 160))
            significance_tab_layout.addWidget(significance_table)
            response_dynamics_tabwidget.addTab(significance_tab, "Significance")

            # Apply the same table style used by residues_conservation_tableWidget
            table_style = ""
            if hasattr(self, 'residues_conservation_tableWidget') and self.residues_conservation_tableWidget is not None:
                table_style = self.residues_conservation_tableWidget.styleSheet()

            if table_style:
                for table_widget in [
                    metrics_table,
                    qc_table,
                    provenance_table,
                    residue_response_table,
                    domain_summary_table,
                    pathway_summary_table,
                    critical_residue_table,
                    network_summary_table,
                    motif_summary_table,
                    significance_table,
                ]:
                    table_widget.setStyleSheet(table_style)

            response_dynamics_tab.setLayout(response_dynamics_tab_layout)

            # --> Dissipation Widget
            dissipation_curve_widget = WidgetPlot(self)

            toolbarSizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
            dissipation_curve_widget.toolbar.setSizePolicy(toolbarSizePolicy)

            dissipationVerticalLayout = QtWidgets.QVBoxLayout()

            dissipationVerticalLayout.addWidget(dissipation_curve_widget.toolbar)
            dissipationVerticalLayout.addWidget(dissipation_curve_widget.canvas)
            dissipation_curve_widget.setLayout(dissipationVerticalLayout)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
            dissipation_curve_widget.setSizePolicy(sizePolicy)
            dissipation_curve_widget.setMinimumSize(QtCore.QSize(0, 400))
            dissipation_curve_widget.setMaximumSize(QtCore.QSize(520, 1200))
            dissipation_curve_widget.setObjectName("dissipation_curve_widget")
            
            # --- START PLOT SUB-TABS ---
            plot_inner_tabwidget = QtWidgets.QTabWidget(plot_tab)
            plot_tab_layout.addWidget(plot_inner_tabwidget)

            # Sub-Tab 1: Response Time Graph
            res_time_tab = QtWidgets.QWidget()
            res_time_layout = QtWidgets.QVBoxLayout(res_time_tab)
            res_time_layout.addWidget(dissipation_curve_widget)
            plot_inner_tabwidget.addTab(res_time_tab, "Response Time Graph")

            # Sub-Tab 2: Propagation Coefficient Plots
            pc_plots_tab = QtWidgets.QWidget()
            pc_plots_layout = QtWidgets.QVBoxLayout(pc_plots_tab)
            plot_inner_tabwidget.addTab(pc_plots_tab, "Propagation Coefficient")
            
            try:
                out_dir = getattr(self, '_network_output_directory', getattr(self, 'output_directory', ''))
                if out_dir and os.path.exists(out_dir):
                    import glob
                    import pandas as pd
                    pc_metric_files = glob.glob(os.path.join(out_dir, "PC_Score_*_propagation_metrics.csv"))

                    if pc_metric_files:
                        stacked_widget = QtWidgets.QStackedWidget()
                        plot_widgets = []
                        plot_titles = []
                        for csv_file in pc_metric_files:
                            try:
                                df = pd.read_csv(csv_file)
                                df_plot = df[df['Propagation_Coefficient (PC)'] > 0].copy()
                                df_plot = df_plot.sort_values('Propagation_Coefficient (PC)', ascending=False).head(30)
                                if df_plot.empty:
                                    continue
                                widget = WidgetPlot()
                                # Response time graph ile aynı ebatlar ve layout
                                toolbarSizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
                                widget.toolbar.setSizePolicy(toolbarSizePolicy)
                                plot_layout = QtWidgets.QVBoxLayout()
                                plot_layout.addWidget(widget.toolbar)
                                plot_layout.addWidget(widget.canvas)
                                widget.setLayout(plot_layout)
                                sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
                                widget.setSizePolicy(sizePolicy)
                                widget.setMinimumSize(QtCore.QSize(0, 500))
                                widget.setMaximumSize(QtCore.QSize(520, 1200))
                                ax = widget.canvas.figure.subplots()
                                bars = ax.barh(df_plot['Residue_ID'], df_plot['Propagation_Coefficient (PC)'], color='salmon', edgecolor='black')
                                ax.tick_params(axis='y', labelsize=8)
                                ax.set_xlabel('Propagation Coefficient (PC)', fontsize=10, fontweight='bold', labelpad=2)
                                ax.set_ylabel('Residue', fontsize=10, fontweight='bold', labelpad=2)
                                ax.set_title('Signal Propagation Capacity', fontsize=12, fontweight='bold')
                                ax.set_ylim(-0.5, len(df_plot['Residue_ID'])-0.5)
                                # Response time graph ile aynı şekilde tight_layout ve margin
                                try:
                                    widget.canvas.figure.tight_layout(rect=(0.15, 0.05, 0.95, 0.95)) #left, bottom, right, top
                                except Exception:
                                    widget.canvas.figure.tight_layout()
                                widget.canvas.draw()
                                stacked_widget.addWidget(widget)
                                plot_widgets.append(widget)
                                plot_titles.append(os.path.basename(csv_file))
                            except Exception as e:
                                err_lbl = QtWidgets.QLabel(f"Error plotting {os.path.basename(csv_file)}: {e}")
                                stacked_widget.addWidget(err_lbl)
                        # Navigation buttons
                        nav_layout = QtWidgets.QHBoxLayout()
                        prev_btn = QtWidgets.QPushButton('Previous')
                        next_btn = QtWidgets.QPushButton('Next')
                        nav_btn_style = """
                        QPushButton {
                            color: white;
                            border: 2px solid rgb(52, 59, 72);
                            border-radius: 5px;
                            background-color: rgb(110, 105, 225);
                            border-width: 1px;
                            padding: 5px 12px;
                            font-size: 11px;
                        }
                        QPushButton:hover {
                            background-color: rgb(22, 200, 244);
                            border: 2px solid rgb(61, 70, 86);
                        }
                        QPushButton:pressed {
                            background-color: rgb(15, 133, 163);
                            border: 2px solid rgb(43, 50, 61);
                        }
                        """
                        prev_btn.setStyleSheet(nav_btn_style)
                        next_btn.setStyleSheet(nav_btn_style)
                        title_lbl = QtWidgets.QLabel()
                        title_lbl.setAlignment(QtCore.Qt.AlignCenter)
                        title_lbl.setStyleSheet("""
                            QLabel {
                                color: #fff;
                                background: transparent;
                                font-size: 12px;
                                font-weight: bold;
                                padding: 4px 12px;
                                border-radius: 4px;
                                letter-spacing: 0.5px;
                            }
                        """)
                        nav_layout.addWidget(prev_btn)
                        nav_layout.addWidget(title_lbl)
                        nav_layout.addWidget(next_btn)
                        def update_title(idx):
                            if 0 <= idx < len(plot_titles):
                                title_lbl.setText(plot_titles[idx])
                            else:
                                title_lbl.setText("")
                        def goto_prev():
                            idx = stacked_widget.currentIndex()
                            if idx > 0:
                                stacked_widget.setCurrentIndex(idx - 1)
                                update_title(idx - 1)
                        def goto_next():
                            idx = stacked_widget.currentIndex()
                            if idx < stacked_widget.count() - 1:
                                stacked_widget.setCurrentIndex(idx + 1)
                                update_title(idx + 1)
                        prev_btn.clicked.connect(goto_prev)
                        next_btn.clicked.connect(goto_next)
                        if plot_widgets:
                            stacked_widget.setCurrentIndex(0)
                            update_title(0)
                        pc_plots_layout.addLayout(nav_layout)
                        # stacked_widget'ın genişliği parent ile birlikte büyüsün
                        stacked_widget.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
                        stacked_widget.setMinimumSize(QtCore.QSize(0, 500))
                        stacked_widget.setMaximumSize(QtCore.QSize(520, 1200))
                        pc_plots_layout.addWidget(stacked_widget, stretch=1)
                    else:
                        lbl = QtWidgets.QLabel("No PC metrics found.")
                        lbl.setAlignment(QtCore.Qt.AlignCenter)
                        pc_plots_layout.addWidget(lbl)
            except Exception as e:
                lbl = QtWidgets.QLabel(f"Error loading PC plots: {e}")
                pc_plots_layout.addWidget(lbl)
            # --- END PLOT SUB-TABS ---

            analysis_data_tabwidget.addTab(paths_tab, "Paths")
            analysis_data_tabwidget.addTab(plot_tab, "Plot")
            tables_tabwidget.addTab(metrics_tab, "Metrics")
            tables_tabwidget.addTab(qc_tab, "QC")
            tables_tabwidget.addTab(provenance_tab, "Provenance")
            tables_tabwidget.addTab(response_dynamics_tab, "Response Dynamics")
            analysis_data_tabwidget.addTab(tables_tab, "Tables")

            last_inner_tab = int(getattr(self, '_analysis_data_tab_last_index', 0))
            if 0 <= last_inner_tab < analysis_data_tabwidget.count():
                analysis_data_tabwidget.setCurrentIndex(last_inner_tab)

            analysis_data_tabwidget.currentChanged.connect(
                lambda idx: setattr(self, '_analysis_data_tab_last_index', idx)
            )

            gridLayout.addWidget(analysis_data_tabwidget, 0, 0, 3, 1)

            ############################################################################################################
            show_navigation_button = QtWidgets.QPushButton(tab)
            show_navigation_button.setEnabled(True)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Preferred)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(show_navigation_button.sizePolicy().hasHeightForWidth())
            show_navigation_button.setSizePolicy(sizePolicy)
            show_navigation_button.setMaximumSize(QtCore.QSize(20, 61))
            show_navigation_button.setStyleSheet(" QPushButton \n"
                                                 "        {\n"
                                                 "         color: white; \n"
                                                 "        border: 2px solid rgb(52, 59, 72); \n"
                                                 "        border-radius: 5px; \n"
                                                 "        background-color:  rgb(255, 17, 100); \n"
                                                 "        border-width: 1px; \n"
                                                 "        outline: none;    \n"
                                                 "        background-color: rgb(110, 105, 225);\n"
                                                 "    \n"
                                                 "        }\n"
                                                 "\n"
                                                 "QPushButton:hover \n"
                                                 "       { \n"
                                                 "        background-color: rgb(22, 200, 244); \n"
                                                 "        border: 2px solid rgb(61, 70, 86);\n"
                                                 "        }\n"
                                                 "\n"
                                                 "QPushButton:pressed \n"
                                                 "        { \n"
                                                 "         background-color:  rgb(15, 133, 163); \n"
                                                 "         border: 2px solid rgb(43, 50, 61);\n"
                                                 "         }      ")
            show_navigation_button.setText("")
            icon11 = QtGui.QIcon()
            icon11.addPixmap(QtGui.QPixmap(":/16x16/icons/16x16/cil-chevron-circle-right-alt.png"), QtGui.QIcon.Normal,
                             QtGui.QIcon.Off)
            show_navigation_button.setIcon(icon11)
            show_navigation_button.setObjectName("show_navigation_button")
            gridLayout.addWidget(show_navigation_button, 0, 1, 3, 1)
            # 0  2  6  1
            hide_navigation_button = QtWidgets.QPushButton(tab)
            hide_navigation_button.setMaximumSize(QtCore.QSize(20, 61))
            hide_navigation_button.setStyleSheet(" QPushButton \n"
                                                 "        {\n"
                                                 "         color: white; \n"
                                                 "        border: 2px solid rgb(52, 59, 72); \n"
                                                 "        border-radius: 5px; \n"
                                                 "        background-color:  rgb(255, 17, 100); \n"
                                                 "        border-width: 1px; \n"
                                                 "        outline: none;    \n"
                                                 "        background-color: rgb(110, 105, 225);\n"
                                                 "    \n"
                                                 "        }\n"
                                                 "\n"
                                                 "QPushButton:hover \n"
                                                 "       { \n"
                                                 "        background-color: rgb(22, 200, 244); \n"
                                                 "        border: 2px solid rgb(61, 70, 86);\n"
                                                 "        }\n"
                                                 "\n"
                                                 "QPushButton:pressed \n"
                                                 "        { \n"
                                                 "         background-color:  rgb(15, 133, 163); \n"
                                                 "         border: 2px solid rgb(43, 50, 61);\n"
                                                 "         }     ")
            hide_navigation_button.setText("")
            hide_navigation_button.setIcon(icon11)
            hide_navigation_button.setObjectName("hide_navigation")
            gridLayout.addWidget(hide_navigation_button, 0, 3, 3, 1)

            analysis_settings_groupBox = QtWidgets.QGroupBox(tab)
            analysis_settings_groupBox.setTitle('Visualization Settings')
            analysis_settings_groupBox.setMinimumSize(QtCore.QSize(170, 0))
            analysis_settings_groupBox.setMaximumSize(QtCore.QSize(170, 16777215))
            analysis_settings_groupBox.setStyleSheet("QGroupBox\n"
                                                     "{\n"
                                                     "    border: 1px solid black;\n"
                                                     "    border-radius: 5px;\n"
                                                     "    border-top-color: rgb(157, 90, 198);\n"
                                                     "    border-left-color: rgb(157, 90, 198);\n"
                                                     "    border-bottom-color: rgb(157, 90, 198);\n"
                                                     "    border-right-color: rgb(157, 90, 198);\n"
                                                     " \n"
                                                     "}")
            analysis_settings_groupBox.setAlignment(
                QtCore.Qt.AlignLeading | QtCore.Qt.AlignLeft | QtCore.Qt.AlignVCenter)
            analysis_settings_groupBox.setObjectName("analysis_settings_groupBox")

            gridLayoutWidget_on_analysis = QtWidgets.QWidget(analysis_settings_groupBox)
            gridLayoutWidget_on_analysis.setGeometry(QtCore.QRect(11, 50, 151, 251))
            gridLayoutWidget_on_analysis.setObjectName("gridLayoutWidget_on_analysis")
            verticalLayout_analysis = QtWidgets.QVBoxLayout(gridLayoutWidget_on_analysis)
            verticalLayout_analysis.setContentsMargins(0, 0, 0, 0)
            verticalLayout_analysis.setObjectName("verticalLayout_analysis")

            activate_pymol_navigation_on_analysis = QtWidgets.QPushButton(gridLayoutWidget_on_analysis)
            activate_pymol_navigation_on_analysis.setText('Activate')
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(activate_pymol_navigation_on_analysis.sizePolicy().hasHeightForWidth())
            activate_pymol_navigation_on_analysis.setSizePolicy(sizePolicy)
            activate_pymol_navigation_on_analysis.setMinimumSize(QtCore.QSize(149, 33))
            activate_pymol_navigation_on_analysis.setMaximumSize(QtCore.QSize(110, 33))
            font = QtGui.QFont()
            font.setFamily("Segoe UI")
            font.setPointSize(9)
            activate_pymol_navigation_on_analysis.setFont(font)
            activate_pymol_navigation_on_analysis.setStyleSheet("            QPushButton \n"
                                                                "            {\n"
                                                                "                color: white; \n"
                                                                "                border: 2px solid rgb(52, 59, 72); \n"
                                                                "                border-radius: 5px; \n"
                                                                "                background-color:  rgb(22, 200, 244); \n"
                                                                "                margin-top:1px; \n"
                                                                "                margin-bottom: 1px; \n"
                                                                "                border-width: 1px; \n"
                                                                "                padding: 5px; \n"
                                                                "                outline: none;\n"
                                                                "            }\n"
                                                                "\n"
                                                                "            QPushButton:hover \n"
                                                                "            { \n"
                                                                "                background-color: rgb(255, 17, 100); \n"
                                                                "                border: 2px solid rgb(61, 70, 86);\n"
                                                                "            }\n"
                                                                "\n"
                                                                "            QPushButton:pressed \n"
                                                                "            { \n"
                                                                "                background-color:  rgb(15, 133, 163); \n"
                                                                "                border: 2px solid rgb(43, 50, 61);\n"
                                                                "            }")
            icon12 = QtGui.QIcon()
            icon12.addPixmap(QtGui.QPixmap(":/24x24/icons/24x24/cil-cursor.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
            activate_pymol_navigation_on_analysis.setIcon(icon12)
            activate_pymol_navigation_on_analysis.setObjectName("activate_pymol_navigation_on_analysis")
            verticalLayout_analysis.addWidget(activate_pymol_navigation_on_analysis)

            deactivate_pymol_navigation_on_analysis = QtWidgets.QPushButton(gridLayoutWidget_on_analysis)
            deactivate_pymol_navigation_on_analysis.setText('Deactivate')
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(deactivate_pymol_navigation_on_analysis.sizePolicy().hasHeightForWidth())
            deactivate_pymol_navigation_on_analysis.setSizePolicy(sizePolicy)
            deactivate_pymol_navigation_on_analysis.setMinimumSize(QtCore.QSize(149, 33))
            deactivate_pymol_navigation_on_analysis.setMaximumSize(QtCore.QSize(110, 33))
            font = QtGui.QFont()
            font.setFamily("Segoe UI")
            font.setPointSize(9)
            deactivate_pymol_navigation_on_analysis.setFont(font)
            deactivate_pymol_navigation_on_analysis.setStyleSheet("            QPushButton \n"
                                                                  "            {\n"
                                                                  "                color: white; \n"
                                                                  "                border: 2px solid rgb(52, 59, 72); \n"
                                                                  "                border-radius: 5px; \n"
                                                                  "                background-color:  rgb(22, 200, 244); \n"
                                                                  "                margin-top:1px; \n"
                                                                  "                margin-bottom: 1px; \n"
                                                                  "                border-width: 1px; \n"
                                                                  "                padding: 5px; \n"
                                                                  "                outline: none;\n"
                                                                  "            }\n"
                                                                  "\n"
                                                                  "            QPushButton:hover \n"
                                                                  "            { \n"
                                                                  "                background-color: rgb(255, 17, 100); \n"
                                                                  "                border: 2px solid rgb(61, 70, 86);\n"
                                                                  "            }\n"
                                                                  "\n"
                                                                  "            QPushButton:pressed \n"
                                                                  "            { \n"
                                                                  "                background-color:  rgb(15, 133, 163); \n"
                                                                  "                border: 2px solid rgb(43, 50, 61);\n"
                                                                  "            }")
            icon13 = QtGui.QIcon()
            icon13.addPixmap(QtGui.QPixmap(":/20x20/icons/20x20/cil-x.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
            deactivate_pymol_navigation_on_analysis.setIcon(icon13)
            deactivate_pymol_navigation_on_analysis.setObjectName("deactivate_pymol_navigation")
            verticalLayout_analysis.addWidget(deactivate_pymol_navigation_on_analysis)

            refresh_pushButton_on_analysis = QtWidgets.QPushButton(gridLayoutWidget_on_analysis)
            refresh_pushButton_on_analysis.setText('Refresh')
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(refresh_pushButton_on_analysis.sizePolicy().hasHeightForWidth())
            refresh_pushButton_on_analysis.setSizePolicy(sizePolicy)
            refresh_pushButton_on_analysis.setMinimumSize(QtCore.QSize(149, 33))
            refresh_pushButton_on_analysis.setMaximumSize(QtCore.QSize(110, 33))
            font = QtGui.QFont()
            font.setFamily("Segoe UI")
            font.setPointSize(9)
            refresh_pushButton_on_analysis.setFont(font)
            refresh_pushButton_on_analysis.setStyleSheet("            QPushButton \n"
                                                         "            {\n"
                                                         "                color: white; \n"
                                                         "                border: 2px solid rgb(52, 59, 72); \n"
                                                         "                border-radius: 5px; \n"
                                                         "                background-color:  rgb(22, 200, 244); \n"
                                                         "                margin-top:1px; \n"
                                                         "                margin-bottom: 1px; \n"
                                                         "                border-width: 1px; \n"
                                                         "                padding: 5px; \n"
                                                         "                outline: none;\n"
                                                         "            }\n"
                                                         "\n"
                                                         "            QPushButton:hover \n"
                                                         "            { \n"
                                                         "                background-color: rgb(255, 17, 100); \n"
                                                         "                border: 2px solid rgb(61, 70, 86);\n"
                                                         "            }\n"
                                                         "\n"
                                                         "            QPushButton:pressed \n"
                                                         "            { \n"
                                                         "                background-color:  rgb(15, 133, 163); \n"
                                                         "                border: 2px solid rgb(43, 50, 61);\n"
                                                         "            }")
            icon14 = QtGui.QIcon()
            icon14.addPixmap(QtGui.QPixmap(":/16x16/icons/16x16/cil-reload.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
            refresh_pushButton_on_analysis.setIcon(icon14)
            refresh_pushButton_on_analysis.setObjectName("refresh_pushButton_on_analysis")
            verticalLayout_analysis.addWidget(refresh_pushButton_on_analysis)

            color_response_panel_pushButton = QtWidgets.QPushButton(gridLayoutWidget_on_analysis)
            color_response_panel_pushButton.setText('Color by Response')
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(color_response_panel_pushButton.sizePolicy().hasHeightForWidth())
            color_response_panel_pushButton.setSizePolicy(sizePolicy)
            color_response_panel_pushButton.setMinimumSize(QtCore.QSize(149, 33))
            color_response_panel_pushButton.setMaximumSize(QtCore.QSize(110, 33))
            font = QtGui.QFont()
            font.setFamily("Segoe UI")
            font.setPointSize(9)
            color_response_panel_pushButton.setFont(font)
            color_response_panel_pushButton.setStyleSheet("            QPushButton \n"
                                                          "            {\n"
                                                          "                color: white; \n"
                                                          "                border: 2px solid rgb(52, 59, 72); \n"
                                                          "                border-radius: 5px; \n"
                                                          "                background-color:  rgb(22, 200, 244); \n"
                                                          "                margin-top:1px; \n"
                                                          "                margin-bottom: 1px; \n"
                                                          "                border-width: 1px; \n"
                                                          "                padding: 5px; \n"
                                                          "                outline: none;\n"
                                                          "            }\n"
                                                          "\n"
                                                          "            QPushButton:hover \n"
                                                          "            { \n"
                                                          "                background-color: rgb(255, 17, 100); \n"
                                                          "                border: 2px solid rgb(61, 70, 86);\n"
                                                          "            }\n"
                                                          "\n"
                                                          "            QPushButton:pressed \n"
                                                          "            { \n"
                                                          "                background-color:  rgb(15, 133, 163); \n"
                                                          "                border: 2px solid rgb(43, 50, 61);\n"
                                                          "            }")
            icon_response = QtGui.QIcon()
            icon_response.addPixmap(QtGui.QPixmap(":/16x16/icons/16x16/cil-brush-alt.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
            color_response_panel_pushButton.setIcon(icon_response)
            color_response_panel_pushButton.setObjectName("color_response_panel_pushButton")
            verticalLayout_analysis.addWidget(color_response_panel_pushButton)

            ss_beatiful_snapshoot_on_analysis = QtWidgets.QPushButton(gridLayoutWidget_on_analysis)
            ss_beatiful_snapshoot_on_analysis.setText('Beautiful Snap')
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(ss_beatiful_snapshoot_on_analysis.sizePolicy().hasHeightForWidth())
            ss_beatiful_snapshoot_on_analysis.setSizePolicy(sizePolicy)
            ss_beatiful_snapshoot_on_analysis.setMinimumSize(QtCore.QSize(149, 33))
            ss_beatiful_snapshoot_on_analysis.setMaximumSize(QtCore.QSize(110, 33))
            font = QtGui.QFont()
            font.setFamily("Segoe UI")
            font.setPointSize(9)
            ss_beatiful_snapshoot_on_analysis.setFont(font)
            ss_beatiful_snapshoot_on_analysis.setStyleSheet("            QPushButton \n"
                                                            "            {\n"
                                                            "                color: white; \n"
                                                            "                border: 2px solid rgb(52, 59, 72); \n"
                                                            "                border-radius: 5px; \n"
                                                            "                background-color:  rgb(22, 200, 244); \n"
                                                            "                margin-top:1px; \n"
                                                            "                margin-bottom: 1px; \n"
                                                            "                border-width: 1px; \n"
                                                            "                padding: 5px; \n"
                                                            "                outline: none;\n"
                                                            "            }\n"
                                                            "\n"
                                                            "            QPushButton:hover \n"
                                                            "            { \n"
                                                            "                background-color: rgb(255, 17, 100); \n"
                                                            "                border: 2px solid rgb(61, 70, 86);\n"
                                                            "            }\n"
                                                            "\n"
                                                            "            QPushButton:pressed \n"
                                                            "            { \n"
                                                            "                background-color:  rgb(15, 133, 163); \n"
                                                            "                border: 2px solid rgb(43, 50, 61);\n"
                                                            "            }")
            icon15 = QtGui.QIcon()
            icon15.addPixmap(QtGui.QPixmap(":/16x16/icons/16x16/cil-camera.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
            ss_beatiful_snapshoot_on_analysis.setIcon(icon15)
            ss_beatiful_snapshoot_on_analysis.setObjectName("ss_beatiful_snapshoot_on_analysis")
            verticalLayout_analysis.addWidget(ss_beatiful_snapshoot_on_analysis)

            save_as_png_on_analysis_pushButton = QtWidgets.QPushButton(gridLayoutWidget_on_analysis)
            save_as_png_on_analysis_pushButton.setText('Save as png')
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(save_as_png_on_analysis_pushButton.sizePolicy().hasHeightForWidth())
            save_as_png_on_analysis_pushButton.setSizePolicy(sizePolicy)
            save_as_png_on_analysis_pushButton.setMinimumSize(QtCore.QSize(149, 33))
            save_as_png_on_analysis_pushButton.setMaximumSize(QtCore.QSize(110, 33))
            font = QtGui.QFont()
            font.setFamily("Segoe UI")
            font.setPointSize(9)
            save_as_png_on_analysis_pushButton.setFont(font)
            save_as_png_on_analysis_pushButton.setStyleSheet("            QPushButton \n"
                                                             "            {\n"
                                                             "                color: white; \n"
                                                             "                border: 2px solid rgb(52, 59, 72); \n"
                                                             "                border-radius: 5px; \n"
                                                             "                background-color:  rgb(22, 200, 244); \n"
                                                             "                margin-top:1px; \n"
                                                             "                margin-bottom: 1px; \n"
                                                             "                border-width: 1px; \n"
                                                             "                padding: 5px; \n"
                                                             "                outline: none;\n"
                                                             "            }\n"
                                                             "\n"
                                                             "            QPushButton:hover \n"
                                                             "            { \n"
                                                             "                background-color: rgb(255, 17, 100); \n"
                                                             "                border: 2px solid rgb(61, 70, 86);\n"
                                                             "            }\n"
                                                             "\n"
                                                             "            QPushButton:pressed \n"
                                                             "            { \n"
                                                             "                background-color:  rgb(15, 133, 163); \n"
                                                             "                border: 2px solid rgb(43, 50, 61);\n"
                                                             "            }")
            icon16 = QtGui.QIcon()
            icon16.addPixmap(QtGui.QPixmap(":/16x16/icons/16x16/cil-data-transfer-down.png"), QtGui.QIcon.Normal,
                             QtGui.QIcon.Off)
            save_as_png_on_analysis_pushButton.setIcon(icon16)
            save_as_png_on_analysis_pushButton.setObjectName("save_as_png_on_analysis_pushButton")
            verticalLayout_analysis.addWidget(save_as_png_on_analysis_pushButton)

            figure_settings_on_analysis_groupBox = QtWidgets.QGroupBox(analysis_settings_groupBox)
            figure_settings_on_analysis_groupBox.setTitle('Figure Settings')
            figure_settings_on_analysis_groupBox.setGeometry(QtCore.QRect(-1, 300, 171, 401))
            figure_settings_on_analysis_groupBox.setObjectName("figure_settings_on_analysis_groupBox")
            verticalLayoutWidget_on_analysis_2 = QtWidgets.QWidget(figure_settings_on_analysis_groupBox)
            verticalLayoutWidget_on_analysis_2.setGeometry(QtCore.QRect(10, 54, 151, 311))
            verticalLayoutWidget_on_analysis_2.setObjectName("verticalLayoutWidget_on_analysis_2")
            figure_settings_on_analysis_verticalLayout = QtWidgets.QVBoxLayout(verticalLayoutWidget_on_analysis_2)
            figure_settings_on_analysis_verticalLayout.setContentsMargins(0, 0, 0, 0)
            figure_settings_on_analysis_verticalLayout.setObjectName("figure_settings_on_analysis_verticalLayout")

            pymol_width_label_on_analysis = QtWidgets.QLabel(verticalLayoutWidget_on_analysis_2)
            pymol_width_label_on_analysis.setText("Width: 1200")
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(pymol_width_label_on_analysis.sizePolicy().hasHeightForWidth())
            pymol_width_label_on_analysis.setSizePolicy(sizePolicy)
            pymol_width_label_on_analysis.setMinimumSize(QtCore.QSize(50, 22))
            pymol_width_label_on_analysis.setMaximumSize(QtCore.QSize(16777215, 22))
            pymol_width_label_on_analysis.setStyleSheet("QLabel {\n"
                                                        "    background-color: rgb(27, 29, 35);\n"
                                                        "    border-radius: 5px;\n"
                                                        "    border: 2px solid rgb(27, 29, 35);\n"
                                                        "    padding: 1px 1px 1px 1px;\n"
                                                        "    border-bottom-color: rgb(157, 90, 198);\n"
                                                        "}\n"
                                                        "\n"
                                                        "\n"
                                                        "QLabel:hover{\n"
                                                        "    border: 2px solid rgb(64, 71, 88);\n"
                                                        "    selection-color: rgb(127, 5, 64);\n"
                                                        "\n"
                                                        "}")
            pymol_width_label_on_analysis.setObjectName("pymol_width_label_on_analysis")
            figure_settings_on_analysis_verticalLayout.addWidget(pymol_width_label_on_analysis)

            width_horizontalSlider_on_analysis = QtWidgets.QSlider(verticalLayoutWidget_on_analysis_2)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(width_horizontalSlider_on_analysis.sizePolicy().hasHeightForWidth())
            width_horizontalSlider_on_analysis.setSizePolicy(sizePolicy)
            width_horizontalSlider_on_analysis.setStyleSheet("\n"
                                                             "\n"
                                                             "QSlider::handle:horizontal \n"
                                                             "    {\n"
                                                             "    background-color:  rgb(255, 17, 100);\n"
                                                             "    border: 2px solid;\n"
                                                             "    width: 8px;\n"
                                                             "    margin: -15px 0px;\n"
                                                             "    }\n"
                                                             "\n"
                                                             "QSlider:horizontal \n"
                                                             "{\n"
                                                             "    min-height: 20px;\n"
                                                             "}\n"
                                                             "\n"
                                                             "QSlider::groove:horizontal \n"
                                                             "{\n"
                                                             "    height: 1px;\n"
                                                             "    background-color: rgb(110, 105, 225);\n"
                                                             "    border: 1px solid;\n"
                                                             "    height: 5px;\n"
                                                             "    margin: 0px;\n"
                                                             "    border-radius: 5px;\n"
                                                             "}\n"
                                                             "\n"
                                                             "QSlider::handle:horizontal \n"
                                                             "{\n"
                                                             "    width: 10px;\n"
                                                             "    margin-top: -10px;\n"
                                                             "    margin-bottom: -10px;\n"
                                                             "    border-radius: 5px;\n"
                                                             "    background-color: rgb(255, 17, 100);\n"
                                                             "    border: 2px solid;\n"
                                                             "}\n"
                                                             "\n"
                                                             "QSlider::handle:horizontal:hover\n"
                                                             "{\n"
                                                             "    background-color: rgb(22, 200, 244);\n"
                                                             "}\n"
                                                             "")
            width_horizontalSlider_on_analysis.setMinimum(800)
            width_horizontalSlider_on_analysis.setMaximum(1920)
            width_horizontalSlider_on_analysis.setSingleStep(20)
            width_horizontalSlider_on_analysis.setPageStep(20)
            width_horizontalSlider_on_analysis.setProperty("value", 1200)
            width_horizontalSlider_on_analysis.setOrientation(QtCore.Qt.Horizontal)
            width_horizontalSlider_on_analysis.setInvertedAppearance(False)
            width_horizontalSlider_on_analysis.setInvertedControls(False)
            width_horizontalSlider_on_analysis.setTickPosition(QtWidgets.QSlider.NoTicks)
            width_horizontalSlider_on_analysis.setTickInterval(20)
            width_horizontalSlider_on_analysis.setObjectName("width_horizontalSlider_on_analysis")
            figure_settings_on_analysis_verticalLayout.addWidget(width_horizontalSlider_on_analysis)

            pymol_height_label_on_analysis = QtWidgets.QLabel(verticalLayoutWidget_on_analysis_2)
            pymol_height_label_on_analysis.setText("Height: 1080")
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(pymol_height_label_on_analysis.sizePolicy().hasHeightForWidth())
            pymol_height_label_on_analysis.setSizePolicy(sizePolicy)
            pymol_height_label_on_analysis.setMinimumSize(QtCore.QSize(50, 22))
            pymol_height_label_on_analysis.setMaximumSize(QtCore.QSize(16777215, 22))
            pymol_height_label_on_analysis.setStyleSheet("QLabel {\n"
                                                         "    background-color: rgb(27, 29, 35);\n"
                                                         "    border-radius: 5px;\n"
                                                         "    border: 2px solid rgb(27, 29, 35);\n"
                                                         "    padding: 1px 1px 1px 1px;\n"
                                                         "    \n"
                                                         "    border-bottom-color: rgb(157, 90, 198);\n"
                                                         "}\n"
                                                         "\n"
                                                         "\n"
                                                         "QLabel:hover{\n"
                                                         "    border: 2px solid rgb(64, 71, 88);\n"
                                                         "    selection-color: rgb(127, 5, 64);\n"
                                                         "\n"
                                                         "}")
            pymol_height_label_on_analysis.setObjectName("pymol_height_label_on_analysis")
            figure_settings_on_analysis_verticalLayout.addWidget(pymol_height_label_on_analysis)

            height_horizontalSlider_on_analysis = QtWidgets.QSlider(verticalLayoutWidget_on_analysis_2)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(height_horizontalSlider_on_analysis.sizePolicy().hasHeightForWidth())
            height_horizontalSlider_on_analysis.setSizePolicy(sizePolicy)
            height_horizontalSlider_on_analysis.setStyleSheet("QSlider::handle:horizontal \n"
                                                              "    {\n"
                                                              "    background-color:  rgb(255, 17, 100);\n"
                                                              "    border: 2px solid;\n"
                                                              "    width: 8px;\n"
                                                              "    margin: -15px 0px;\n"
                                                              "    }\n"
                                                              "\n"
                                                              "QSlider:horizontal \n"
                                                              "{\n"
                                                              "    min-height: 20px;\n"
                                                              "}\n"
                                                              "\n"
                                                              "QSlider::groove:horizontal \n"
                                                              "{\n"
                                                              "    height: 1px;\n"
                                                              "    background-color: rgb(110, 105, 225);\n"
                                                              "    border: 1px solid;\n"
                                                              "    height: 5px;\n"
                                                              "    margin: 0px;\n"
                                                              "    border-radius: 5px;\n"
                                                              "}\n"
                                                              "\n"
                                                              "QSlider::handle:horizontal \n"
                                                              "{\n"
                                                              "    width: 10px;\n"
                                                              "    margin-top: -10px;\n"
                                                              "    margin-bottom: -10px;\n"
                                                              "    border-radius: 5px;\n"
                                                              "    background-color: rgb(255, 17, 100);\n"
                                                              "    border: 2px solid;\n"
                                                              "}\n"
                                                              "\n"
                                                              "QSlider::handle:horizontal:hover\n"
                                                              "{\n"
                                                              "    background-color: rgb(22, 200, 244);\n"
                                                              "}\n"
                                                              "")
            height_horizontalSlider_on_analysis.setMinimum(600)
            height_horizontalSlider_on_analysis.setMaximum(1080)
            height_horizontalSlider_on_analysis.setSingleStep(20)
            height_horizontalSlider_on_analysis.setPageStep(20)
            height_horizontalSlider_on_analysis.setProperty("value", 1080)
            height_horizontalSlider_on_analysis.setOrientation(QtCore.Qt.Horizontal)
            height_horizontalSlider_on_analysis.setTickInterval(20)
            height_horizontalSlider_on_analysis.setObjectName("height_horizontalSlider_on_analysis")
            figure_settings_on_analysis_verticalLayout.addWidget(height_horizontalSlider_on_analysis)

            pymol_dpi_label_on_analysis = QtWidgets.QLabel(verticalLayoutWidget_on_analysis_2)
            pymol_dpi_label_on_analysis.setText("Dpi: 150")
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(pymol_dpi_label_on_analysis.sizePolicy().hasHeightForWidth())
            pymol_dpi_label_on_analysis.setSizePolicy(sizePolicy)
            pymol_dpi_label_on_analysis.setMinimumSize(QtCore.QSize(50, 22))
            pymol_dpi_label_on_analysis.setMaximumSize(QtCore.QSize(16777215, 22))
            pymol_dpi_label_on_analysis.setStyleSheet("QLabel {\n"
                                                      "    background-color: rgb(27, 29, 35);\n"
                                                      "    border-radius: 5px;\n"
                                                      "    border: 2px solid rgb(27, 29, 35);\n"
                                                      "    padding: 1px 1px 1px 1px;\n"
                                                      "    \n"
                                                      "    border-bottom-color: rgb(157, 90, 198);\n"
                                                      "}\n"
                                                      "\n"
                                                      "\n"
                                                      "QLabel:hover{\n"
                                                      "    border: 2px solid rgb(64, 71, 88);\n"
                                                      "    selection-color: rgb(127, 5, 64);\n"
                                                      "\n"
                                                      "}")
            pymol_dpi_label_on_analysis.setObjectName("pymol_dpi_label_on_analysis")
            figure_settings_on_analysis_verticalLayout.addWidget(pymol_dpi_label_on_analysis)

            dpi_horizontalSlider_on_analysis = QtWidgets.QSlider(verticalLayoutWidget_on_analysis_2)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(dpi_horizontalSlider_on_analysis.sizePolicy().hasHeightForWidth())
            dpi_horizontalSlider_on_analysis.setSizePolicy(sizePolicy)
            dpi_horizontalSlider_on_analysis.setStyleSheet("QSlider::handle:horizontal \n"
                                                           "    {\n"
                                                           "    background-color:  rgb(255, 17, 100);\n"
                                                           "    border: 2px solid;\n"
                                                           "    width: 8px;\n"
                                                           "    margin: -15px 0px;\n"
                                                           "    }\n"
                                                           "\n"
                                                           "QSlider:horizontal \n"
                                                           "{\n"
                                                           "    min-height: 20px;\n"
                                                           "}\n"
                                                           "\n"
                                                           "QSlider::groove:horizontal \n"
                                                           "{\n"
                                                           "    height: 1px;\n"
                                                           "    background-color: rgb(110, 105, 225);\n"
                                                           "    border: 1px solid;\n"
                                                           "    height: 5px;\n"
                                                           "    margin: 0px;\n"
                                                           "    border-radius: 5px;\n"
                                                           "}\n"
                                                           "\n"
                                                           "QSlider::handle:horizontal \n"
                                                           "{\n"
                                                           "    width: 10px;\n"
                                                           "    margin-top: -10px;\n"
                                                           "    margin-bottom: -10px;\n"
                                                           "    border-radius: 5px;\n"
                                                           "    background-color: rgb(255, 17, 100);\n"
                                                           "    border: 2px solid;\n"
                                                           "}\n"
                                                           "\n"
                                                           "QSlider::handle:horizontal:hover\n"
                                                           "{\n"
                                                           "    background-color: rgb(22, 200, 244);\n"
                                                           "}\n"
                                                           "")
            dpi_horizontalSlider_on_analysis.setMinimum(100)
            dpi_horizontalSlider_on_analysis.setMaximum(300)
            dpi_horizontalSlider_on_analysis.setSingleStep(50)
            dpi_horizontalSlider_on_analysis.setPageStep(50)
            dpi_horizontalSlider_on_analysis.setProperty("value", 150)
            dpi_horizontalSlider_on_analysis.setOrientation(QtCore.Qt.Horizontal)
            dpi_horizontalSlider_on_analysis.setTickInterval(50)
            dpi_horizontalSlider_on_analysis.setObjectName("dpi_horizontalSlider_on_analysis")
            figure_settings_on_analysis_verticalLayout.addWidget(dpi_horizontalSlider_on_analysis)

            pymol_ray_label_on_analysis = QtWidgets.QLabel(verticalLayoutWidget_on_analysis_2)
            pymol_ray_label_on_analysis.setText("Ray: 1")
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(pymol_ray_label_on_analysis.sizePolicy().hasHeightForWidth())
            pymol_ray_label_on_analysis.setSizePolicy(sizePolicy)
            pymol_ray_label_on_analysis.setMinimumSize(QtCore.QSize(0, 22))
            pymol_ray_label_on_analysis.setMaximumSize(QtCore.QSize(16777215, 22))
            pymol_ray_label_on_analysis.setStyleSheet("QLabel {\n"
                                                      "    background-color: rgb(27, 29, 35);\n"
                                                      "    border-radius: 5px;\n"
                                                      "    border: 2px solid rgb(27, 29, 35);\n"
                                                      "    padding: 1px 1px 1px 1px;\n"
                                                      "    \n"
                                                      "    border-bottom-color: rgb(157, 90, 198);\n"
                                                      "}\n"
                                                      "\n"
                                                      "\n"
                                                      "QLabel:hover{\n"
                                                      "    border: 2px solid rgb(64, 71, 88);\n"
                                                      "    selection-color: rgb(127, 5, 64);\n"
                                                      "\n"
                                                      "}")
            pymol_ray_label_on_analysis.setObjectName("pymol_ray_label_on_analysis")
            figure_settings_on_analysis_verticalLayout.addWidget(pymol_ray_label_on_analysis)

            ray_horizontalSlider_on_analysis = QtWidgets.QSlider(verticalLayoutWidget_on_analysis_2)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(ray_horizontalSlider_on_analysis.sizePolicy().hasHeightForWidth())
            ray_horizontalSlider_on_analysis.setSizePolicy(sizePolicy)
            ray_horizontalSlider_on_analysis.setStyleSheet("QSlider::handle:horizontal \n"
                                                           "    {\n"
                                                           "    background-color:  rgb(255, 17, 100);\n"
                                                           "    border: 2px solid;\n"
                                                           "    width: 8px;\n"
                                                           "    margin: -15px 0px;\n"
                                                           "    }\n"
                                                           "\n"
                                                           "QSlider:horizontal \n"
                                                           "{\n"
                                                           "    min-height: 20px;\n"
                                                           "}\n"
                                                           "\n"
                                                           "QSlider::groove:horizontal \n"
                                                           "{\n"
                                                           "    height: 1px;\n"
                                                           "    background-color: rgb(110, 105, 225);\n"
                                                           "    border: 1px solid;\n"
                                                           "    height: 5px;\n"
                                                           "    margin: 0px;\n"
                                                           "    border-radius: 5px;\n"
                                                           "}\n"
                                                           "\n"
                                                           "QSlider::handle:horizontal \n"
                                                           "{\n"
                                                           "    width: 10px;\n"
                                                           "    margin-top: -10px;\n"
                                                           "    margin-bottom: -10px;\n"
                                                           "    border-radius: 5px;\n"
                                                           "    background-color: rgb(255, 17, 100);\n"
                                                           "    border: 2px solid;\n"
                                                           "}\n"
                                                           "\n"
                                                           "QSlider::handle:horizontal:hover\n"
                                                           "{\n"
                                                           "    background-color: rgb(22, 200, 244);\n"
                                                           "}\n"
                                                           "")
            ray_horizontalSlider_on_analysis.setMinimum(1)
            ray_horizontalSlider_on_analysis.setMaximum(5)
            ray_horizontalSlider_on_analysis.setPageStep(1)
            ray_horizontalSlider_on_analysis.setProperty("value", 1)
            ray_horizontalSlider_on_analysis.setOrientation(QtCore.Qt.Horizontal)
            ray_horizontalSlider_on_analysis.setTickInterval(1)
            ray_horizontalSlider_on_analysis.setObjectName("ray_horizontalSlider_on_analysis")
            figure_settings_on_analysis_verticalLayout.addWidget(ray_horizontalSlider_on_analysis)

            get_figure_pushButton_on_analysis = QtWidgets.QPushButton(verticalLayoutWidget_on_analysis_2)
            get_figure_pushButton_on_analysis.setText("Get Figure")
            get_figure_pushButton_on_analysis.setMinimumSize(QtCore.QSize(0, 33))
            font = QtGui.QFont()
            font.setFamily("Segoe UI")
            font.setPointSize(9)
            get_figure_pushButton_on_analysis.setFont(font)
            get_figure_pushButton_on_analysis.setStyleSheet("QPushButton \n"
                                                            "{\n"
                                                            "    border: 2px solid rgb(52, 59, 72);\n"
                                                            "    border-radius: 5px;        \n"
                                                            "    background-color: rgb(202, 116, 220);\n"
                                                            "    \n"
                                                            "    background-color: rgb(253, 1, 136);\n"
                                                            "}\n"
                                                            "\n"
                                                            "QPushButton:hover \n"
                                                            "{\n"
                                                            "    \n"
                                                            "    background-color: rgb(22, 200, 244); \n"
                                                            "    border: 2px solid rgb(61, 70, 86);\n"
                                                            "}\n"
                                                            "QPushButton:pressed \n"
                                                            "{    \n"
                                                            "    background-color: rgb(255, 170, 0);\n"
                                                            "    border: 2px solid rgb(43, 50, 61);\n"
                                                            "}")
            icon17 = QtGui.QIcon()
            icon17.addPixmap(QtGui.QPixmap(":/16x16/icons/16x16/cil-image-plus.png"), QtGui.QIcon.Normal,
                             QtGui.QIcon.Off)
            get_figure_pushButton_on_analysis.setIcon(icon17)
            get_figure_pushButton_on_analysis.setObjectName("get_figure_pushButton_on_analysis")
            figure_settings_on_analysis_verticalLayout.addWidget(get_figure_pushButton_on_analysis)

            hide_figure_settings_on_analysis_pushButton = QtWidgets.QPushButton(figure_settings_on_analysis_groupBox)
            hide_figure_settings_on_analysis_pushButton.setEnabled(True)
            hide_figure_settings_on_analysis_pushButton.setGeometry(QtCore.QRect(55, 370, 65, 20))
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(hide_figure_settings_on_analysis_pushButton.sizePolicy().hasHeightForWidth())
            hide_figure_settings_on_analysis_pushButton.setSizePolicy(sizePolicy)
            hide_figure_settings_on_analysis_pushButton.setMinimumSize(QtCore.QSize(65, 20))
            hide_figure_settings_on_analysis_pushButton.setMaximumSize(QtCore.QSize(65, 20))
            font = QtGui.QFont()
            font.setWeight(50)
            font.setBold(False)
            hide_figure_settings_on_analysis_pushButton.setFont(font)
            hide_figure_settings_on_analysis_pushButton.setStyleSheet(" QPushButton \n"
                                                                      "        {\n"
                                                                      "         color: white; \n"
                                                                      "        border: 2px solid rgb(52, 59, 72); \n"
                                                                      "        border-radius: 5px; \n"
                                                                      "        background-color:  rgb(255, 17, 100); \n"
                                                                      "        border-width: 1px; \n"
                                                                      "        outline: none;    \n"
                                                                      "        background-color: rgb(110, 105, 225);\n"
                                                                      "        }\n"
                                                                      "\n"
                                                                      "QPushButton:hover \n"
                                                                      "       { \n"
                                                                      "        background-color: rgb(22, 200, 244); \n"
                                                                      "        border: 2px solid rgb(61, 70, 86);\n"
                                                                      "        }\n"
                                                                      "\n"
                                                                      "QPushButton:pressed \n"
                                                                      "        { \n"
                                                                      "         background-color:  rgb(15, 133, 163); \n"
                                                                      "         border: 2px solid rgb(43, 50, 61);\n"
                                                                      "         }     ")
            hide_figure_settings_on_analysis_pushButton.setText("")
            icon18 = QtGui.QIcon()
            icon18.addPixmap(QtGui.QPixmap(":/16x16/icons/16x16/cil-chevron-circle-up-alt.png"),
                             QtGui.QIcon.Normal, QtGui.QIcon.Off)
            hide_figure_settings_on_analysis_pushButton.setIcon(icon18)
            hide_figure_settings_on_analysis_pushButton.setFlat(False)
            hide_figure_settings_on_analysis_pushButton.setObjectName("hide_figure_settings_pushButton")

            def _resolve_response_energy_files(response_file_path):
                normalized_response_path = os.path.abspath(os.path.expanduser(response_file_path.strip()))
                response_dir = os.path.dirname(normalized_response_path)
                response_name = os.path.basename(normalized_response_path)

                suffix = ""
                if response_name.startswith('responseTimes_') and response_name.lower().endswith('.csv'):
                    suffix = response_name[len('responseTimes_'):-4]

                reference_candidates = []
                if suffix:
                    reference_candidates.append(os.path.join(response_dir, f'reference_energy_file_{suffix}.csv'))
                reference_candidates.append(os.path.join(response_dir, 'reference_energy_file.csv'))

                modified_candidates = []
                if suffix:
                    modified_candidates.append(os.path.join(response_dir, f'modified_energy_file_{suffix}.csv'))
                modified_candidates.append(os.path.join(response_dir, 'modified_energy_file.csv'))

                discovered_reference = []
                discovered_modified = []
                try:
                    if os.path.isdir(response_dir):
                        all_files = os.listdir(response_dir)
                        discovered_reference = sorted(
                            os.path.join(response_dir, file_name)
                            for file_name in all_files
                            if file_name.lower().startswith('reference_energy_file') and file_name.lower().endswith('.csv')
                        )
                        discovered_modified = sorted(
                            os.path.join(response_dir, file_name)
                            for file_name in all_files
                            if file_name.lower().startswith('modified_energy_file') and file_name.lower().endswith('.csv')
                        )
                except Exception as list_error:
                    import sys
                    print(f"Warning: Could not list directory {response_dir}: {list_error}", file=sys.stderr)

                reference_energy_file = next((path for path in reference_candidates if os.path.exists(path)), None)
                if reference_energy_file is None and discovered_reference:
                    reference_energy_file = discovered_reference[0]

                modified_energy_file = next((path for path in modified_candidates if os.path.exists(path)), None)
                if modified_energy_file is None and discovered_modified:
                    if suffix:
                        suffix_matches = [
                            path for path in discovered_modified
                            if os.path.basename(path).lower() == f'modified_energy_file_{suffix}'.lower() + '.csv'
                        ]
                        if suffix_matches:
                            modified_energy_file = suffix_matches[0]
                    if modified_energy_file is None:
                        modified_energy_file = discovered_modified[0]

                return reference_energy_file, modified_energy_file, discovered_reference, discovered_modified

            def _refresh_response_dynamics_views(response_file_path):
                source_residue_text = self.source_res_comboBox.currentText()
                row, col, response_count, plot_name, fit_curve, metrics = getResponseTimeGraph(response_file_path)

                residue_names = None
                try:
                    pdb_path = self.boundForm_pdb_lineedit.text().strip()
                    if pdb_path and os.path.exists(pdb_path):
                        _, residue_names = get_residues(pdb_path)
                except Exception as residue_name_error:
                    import sys
                    print(f"Warning: Could not derive residue names from topology file: {residue_name_error}", file=sys.stderr)

                try:
                    response_payload = build_response_dynamics_payload(
                        response_file_path,
                        metrics,
                        frame_time_delta=1.0,
                        residue_names=residue_names,
                    )
                    populate_residue_response_table(residue_response_table, response_payload['residue_rows'])
                    populate_domain_summary_table(domain_summary_table, response_payload['domain_rows'])
                    populate_metrics_table(metrics_table, response_payload['metrics_rows'])
                    populate_provenance_table(provenance_table, response_payload['provenance_rows'])
                except Exception as analyzer_error:
                    import sys
                    print(f"Warning: Could not build response dynamics payload: {analyzer_error}", file=sys.stderr)

                if source_residue_text == '':
                    dissipation_curve_widget.canvas.plot(
                        response_count,
                        source_residue=None,
                        plot_name=plot_name,
                        fitted_data=fit_curve,
                    )
                else:
                    dissipation_curve_widget.canvas.plot(
                        response_count,
                        source_residue=source_residue_text,
                        plot_name=plot_name,
                        fitted_data=fit_curve,
                    )

            def _recalculate_response_time_on_analysis():
                selected_response_path = str(self.response_time_lineEdit.text()).strip()
                if not os.path.isabs(selected_response_path):
                    selected_response_path = os.path.abspath(os.path.join(_safe_getcwd(), selected_response_path))
                selected_response_path = os.path.abspath(os.path.expanduser(selected_response_path))

                response_parent_dir = os.path.dirname(selected_response_path)
                if response_parent_dir and not os.path.exists(response_parent_dir):
                    os.makedirs(response_parent_dir, exist_ok=True)

                self.response_time_lineEdit.setText(selected_response_path)

                if not (os.path.exists(selected_response_path) and selected_response_path.lower().endswith('.csv')):
                    Message_Boxes.Warning_message(
                        self,
                        "Response Time File Missing",
                        "Please select a valid response time CSV file before recalculation.",
                        Style.MessageBox_stylesheet,
                    )
                    return

                reference_energy_file, modified_energy_file, discovered_reference, discovered_modified = _resolve_response_energy_files(selected_response_path)
                if reference_energy_file is None or modified_energy_file is None:
                    discovered_reference_text = ', '.join(os.path.basename(path) for path in discovered_reference) or 'none'
                    discovered_modified_text = ', '.join(os.path.basename(path) for path in discovered_modified) or 'none'
                    Message_Boxes.Warning_message(
                        self,
                        "Energy File Missing",
                        "Reference or modified energy file was not found for recalculation.\n"
                        f"Found reference files: {discovered_reference_text}\n"
                        f"Found modified files: {discovered_modified_text}",
                        Style.MessageBox_stylesheet,
                    )
                    return

                response_threshold_value = float(response_threshold_spinbox.value())
                self._analysis_response_threshold = response_threshold_value

                try:
                    try:
                        from no_gui.response_time_creator import get_residue_response_times
                    except ModuleNotFoundError:
                        from mdpertool.no_gui.response_time_creator import get_residue_response_times

                    get_residue_response_times(
                        reference_energy_file,
                        modified_energy_file,
                        output_name=selected_response_path,
                        response_threshold=response_threshold_value,
                    )
                    _refresh_response_dynamics_views(selected_response_path)
                    _refresh_pathway_and_critical_views(show_done_message=False)

                    Message_Boxes.Information_message(
                        self,
                        "Response Recalculated",
                        f"Response time has been recalculated with threshold {response_threshold_value:.6f}.",
                        Style.MessageBox_stylesheet,
                    )
                except FileNotFoundError as recalc_file_error:
                    error_str = str(recalc_file_error)

                    # Try to extract filename from exception
                    missing_path = error_str
                    if "No such file or directory:" in error_str:
                        missing_path = error_str.split("No such file or directory:")[1].strip().strip("'\"")
                    elif "[Errno 2]" in error_str:
                        parts = error_str.split("[Errno 2]")
                        if len(parts) > 1:
                            missing_path = parts[1].strip().strip("'\"")
                    
                    if hasattr(recalc_file_error, 'filename') and recalc_file_error.filename:
                        missing_path = recalc_file_error.filename
                    
                    Message_Boxes.Warning_message(
                        self,
                        "Recalculation Failed",
                        "Response recalculation could not start because a required input file is missing.\n"
                        f"Missing file: {missing_path}\n\n"
                        "Please ensure that both reference_energy_file*.csv and modified_energy_file*.csv exist in the analysis directory.",
                        Style.MessageBox_stylesheet,
                    )
                except Exception as recalc_error:
                    Message_Boxes.Warning_message(
                        self,
                        "Recalculation Failed",
                        str(recalc_error),
                        Style.MessageBox_stylesheet,
                    )

            recalculate_response_button.clicked.connect(_recalculate_response_time_on_analysis)

            possible_path = str(self.response_time_lineEdit.text())
            if os.path.exists(possible_path.strip()) and possible_path.split('.')[-1] == 'csv':
                _refresh_response_dynamics_views(possible_path)

            # ################################# ==> START - 3D WIDGETS LOCATING <== ################################## #
            pyMOL_3D_analysis_frame = QtWidgets.QFrame(tab)
            pyMOL_3D_analysis_frame.setStyleSheet("QFrame {\n"
                                                  "   border: 1px solid black;\n"
                                                  "   border-radius: 5px;\n"
                                                  "   border-top-color: rgb(157, 90, 198);\n"
                                                  "   border-left-color: rgb(157, 90, 198);\n"
                                                  "   border-bottom-color: rgb(157, 90, 198);\n"
                                                  "   border-right-color: rgb(157, 90, 198);\n"
                                                  "   margin-top: 5px;\n"
                                                  "}")

            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(pyMOL_3D_analysis_frame.sizePolicy().hasHeightForWidth())
            pyMOL_3D_analysis_frame.setSizePolicy(sizePolicy)
            pyMOL_3D_analysis_frame.setFrameShape(QtWidgets.QFrame.StyledPanel)
            pyMOL_3D_analysis_frame.setFrameShadow(QtWidgets.QFrame.Raised)
            pyMOL_3D_analysis_frame.setObjectName("pyMOL_3D_analysis_frame")
            gridLayout.addWidget(pyMOL_3D_analysis_frame, 0, 2, 3, 1)
            horizontalLayout.addLayout(gridLayout)
            self.analysis_TabWidget.addTab(tab, "Analysis " + str(self.tab_count_on_analysis))
            horizontalLayout.addWidget(analysis_settings_groupBox)

            Protein3DNetworkView = PymolQtWidget(self)
            Protein3DNetworkView.change_default_background()
            verticalLayoutProteinNetworkView = QVBoxLayout(pyMOL_3D_analysis_frame)
            verticalLayoutProteinNetworkView.addWidget(Protein3DNetworkView)
            self.setLayout(verticalLayoutProteinNetworkView)
            Protein3DNetworkView.loadMolFile(self.boundForm_pdb_lineedit.text())
            Protein3DNetworkView.update()
            Protein3DNetworkView.show()
            verticalLayoutProteinNetworkView.setContentsMargins(0, 0, 0, 0)

            def _apply_response_coloring():
                selected_response_path = str(self.response_time_lineEdit.text()).strip()
                if os.path.exists(selected_response_path) and selected_response_path.lower().endswith('.csv'):
                    Protein3DNetworkView.show_energy_dissipation(response_time_file_path=selected_response_path)
                    return

                Message_Boxes.Warning_message(
                    self,
                    "Response Time File Missing",
                    "Please select a valid response time CSV file before applying 3D coloring.",
                    Style.MessageBox_stylesheet,
                )

            color_response_panel_pushButton.clicked.connect(_apply_response_coloring)
            # ################################# ==> END - 3D WIDGETS LOCATING <== ################################## #
            # matplotlib_widget = WidgetPlot(self)
            # verticalLayout.addWidget(self.matplotlib_widget.toolbar)
            # verticalLayout.addWidget(self.matplotlib_widget.canvas)
            ##########################################################################################################
            Helper_Functions.visualization_Handle_buttons_changing_on_analysis(self, analysis_settings_groupBox,
                                                                               show_navigation_button,
                                                                               hide_navigation_button)
            Helper_Functions.Handle_Buttons_on_analysis(self, analysis_settings_groupBox, show_navigation_button,
                                                        hide_navigation_button)

            activate_pymol_navigation_on_analysis.clicked.connect(
                lambda: Helper_Functions.activate_navigation_on_Pymol(self, Protein3DNetworkView))
            deactivate_pymol_navigation_on_analysis.clicked.connect(
                lambda: Helper_Functions.deactivate_navigation_on_Pymol(self, Protein3DNetworkView))
            refresh_pushButton_on_analysis.clicked.connect(
                lambda: Helper_Functions.clear_residue_labels(self, Protein3DNetworkView))
            ss_beatiful_snapshoot_on_analysis.clicked.connect(
                lambda: Helper_Functions.show_beautiful_in_Pymol(self, Protein3DNetworkView))
            get_figure_pushButton_on_analysis.clicked.connect(
                lambda: Helper_Functions.save_as_png_Pymol(self, Protein3DNetworkView,
                                                           width_horizontalSlider_on_analysis,
                                                           height_horizontalSlider_on_analysis,
                                                           dpi_horizontalSlider_on_analysis,
                                                           ray_horizontalSlider_on_analysis))

            Helper_Functions.Handle_Save_Figure_Options_on_analysis(self, save_as_png_on_analysis_pushButton,
                                                                    hide_figure_settings_on_analysis_pushButton,
                                                                    width_horizontalSlider_on_analysis,
                                                                    height_horizontalSlider_on_analysis,
                                                                    dpi_horizontalSlider_on_analysis,
                                                                    ray_horizontalSlider_on_analysis,
                                                                    figure_settings_on_analysis_groupBox,
                                                                    pymol_width_label_on_analysis,
                                                                    pymol_height_label_on_analysis,
                                                                    pymol_dpi_label_on_analysis,
                                                                    pymol_ray_label_on_analysis)

            Helper_Functions.Handle_Save_Figure_Options_on_analysis_Changed(self, figure_settings_on_analysis_groupBox)
            ##########################################################################################################

            done_message_shown = False
            current_shortest_path_graphs = []

            if self.initial_network is not None and len(self.initial_network.nodes()) > 0:
                global_betweenness = nx.betweenness_centrality(self.initial_network)
            else:
                global_betweenness = {}

            target_res_list, target_graph_list = extract_target_graph_pairs(
                clean_log_list,
                all_graph_list,
                amino_acid_residues,
            )

            def _refresh_pathway_and_critical_views(show_done_message=False):
                nonlocal done_message_shown
                local_source_res = self.source_res_comboBox.currentText().strip()
                pathway_rows_local, residue_path_hits_local, shortest_path_strings_local, all_paths_messages_local = summarize_target_pathways(
                    local_source_res,
                    target_res_list,
                    target_graph_list,
                )

                shortest_path_listWidget.clear()
                current_shortest_path_graphs.clear()
                current_shortest_path_graphs.extend(target_graph_list[:len(shortest_path_strings_local)])

                for cnt, shortest_path_text in enumerate(shortest_path_strings_local):
                    item = QListWidgetItem(shortest_path_text)
                    item.setBackground(QColor(colors[cnt % len(colors)]))
                    shortest_path_listWidget.addItem(item)

                pathway_rows_sorted = sorted(pathway_rows_local, key=lambda row: (row[3] == 'No Path', str(row[2]), row[0]))
                populate_pathway_summary_table(pathway_summary_table, pathway_rows_sorted)
                reachable_count, unreachable_count = count_reachability(pathway_rows_sorted)
                populate_reachability_qc(qc_table, reachable_count, unreachable_count)

                critical_rows_local = build_critical_residue_rows(residue_path_hits_local, global_betweenness, top_n=20)
                populate_critical_residue_table(critical_residue_table, critical_rows_local)

                if show_done_message:
                    all_path_string = build_done_message(list(all_paths_messages_local))
                    Message_Boxes.Information_message(self, "DONE !", all_path_string, Style.MessageBox_stylesheet)
                    done_message_shown = True

            _refresh_pathway_and_critical_views(show_done_message=True)

            try:
                network_summary_rows = build_network_summary_rows(self.initial_network, all_graph_list)
                populate_network_summary_table(network_summary_table, network_summary_rows)
            except Exception as summary_error:
                import sys
                print(f"Warning: Could not build network summary table: {summary_error}", file=sys.stderr)

            try:
                intersection_graph = build_intersection_graph(all_graph_list)
                if intersection_graph.number_of_nodes() > 0:
                    motif_scope = "Intersection"
                    motif_graph = intersection_graph
                else:
                    motif_scope = "Union"
                    motif_graph = build_union_graph(all_graph_list)

                motif_rows = build_motif_summary_rows(
                    graph=motif_graph,
                    scope_name=motif_scope,
                    motif_sizes=(3, 4),
                    max_combinations=400000,
                    top_k_per_size=10,
                )
                populate_motif_summary_table(motif_summary_table, motif_rows)
            except Exception as motif_error:
                import sys
                print(f"Warning: Could not build motif summary table: {motif_error}", file=sys.stderr)

            try:
                significance_graph = intersection_graph if 'intersection_graph' in locals() and intersection_graph.number_of_nodes() > 0 else build_union_graph(all_graph_list)
                significance_rows = build_network_significance_rows(
                    initial_graph=self.initial_network,
                    core_graph=significance_graph,
                    betweenness_rank_size=20,
                )
                populate_significance_table(significance_table, significance_rows)
            except Exception as significance_error:
                import sys
                print(f"Warning: Could not build significance table: {significance_error}", file=sys.stderr)

            def _show_selected_shortest_path(item):
                selected_row = shortest_path_listWidget.currentRow()
                if not (0 <= selected_row < len(current_shortest_path_graphs)):
                    return
                Functions.show_shortest_paths_on_3D_ProteinView(
                    self,
                    item,
                    Protein3DNetworkView,
                    current_shortest_path_graphs[selected_row],
                )

            shortest_path_listWidget.itemDoubleClicked.connect(
                _show_selected_shortest_path)
            # #########################################################################################################

            # #########################################################################################################

            try:
                if len(self.initial_network.nodes()) > 0:
                    try:
                        arrows_cordinates, intersection_node_list = Pymol_Visualize_Path(graph=self.initial_network,
                                                                                         pdb_file=self.pdb)

                        # ------------------> 3D NETWORK VISUALIZATION USING PYMOL / START <---------------------- #
                        Protein3DNetworkView.show_energy_dissipation(response_time_file_path=self.retime_file)

                        for arrow_coord in arrows_cordinates:
                            Protein3DNetworkView.create_interacting_Residues(atom1=arrow_coord[0],
                                                                             atom2=arrow_coord[1],
                                                                             radius=0.05, gap=0.4, hradius=0.4,
                                                                             hlength=0.8, color='cyan')

                        # MAKE PYMOL VISUALIZATION BETTER
                        Protein3DNetworkView._pymol.cmd.set('cartoon_oval_length', 0.8)  # default is 1.20)
                        Protein3DNetworkView._pymol.cmd.set('cartoon_oval_width', 0.2)
                        Protein3DNetworkView._pymol.cmd.center(selection="all", state=0, origin=1, animate=0)
                        Protein3DNetworkView._pymol.cmd.zoom('all', buffer=0.0, state=0, complete=0)
                        Protein3DNetworkView.update()
                        Protein3DNetworkView.show()

                        # --------------------> 2D NETWORK VISUALIZATION USING visJS / START <-------------------- #
                    except Exception as error:
                        import traceback
                        exc_type, exc_obj, exc_tb = sys.exc_info()
                        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
                        lineno = exc_tb.tb_lineno  # Satır numarasını almak için
                        print(f"Error: {error}")
                        print(f"Exception Type: {exc_type}")
                        print(f"File: {fname}")
                        print(f"Line Number: {lineno}")
                        print(f"Traceback: {traceback.format_exc()}")

                    all = ''
                    for j in clean_log_list:
                        # if j != '':
                        all = all + '\n' + j
                    if not done_message_shown:
                        Message_Boxes.Information_message(self, "DONE !", all, Style.MessageBox_stylesheet)
                    del self.log_holder, self.network_holder

                else:
                    Message_Boxes.Information_message(self, "DONE !", "There is no Intersection Network :(",
                                                      Style.MessageBox_stylesheet)
            except Exception as e:
                exc_type, exc_obj, exc_tb = sys.exc_info()
                fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
                print(exc_type, fname, exc_tb.tb_lineno)

        else:
            Message_Boxes.Warning_message(
                self,
                "No Suitable Graph Found",
                "No suitable graph was found for the selected search parameters.\n"
                "Try lowering the response threshold, adjusting cutoff values, or selecting a different source residue.",
                Style.MessageBox_stylesheet,
            )

    def show_shortest_paths_on_3D_ProteinView(self, item, PyMOL_Widget, selected_graph):
        processed_path = [x.strip() for x in item.text().split('-->')]
        shortest_path_arrow_coords = Shortest_Path_Visualize(pdb_file=self.pdb, selected_path=processed_path)
        arrows_cordinates, _ = Pymol_Visualize_Path(graph=selected_graph, pdb_file=self.pdb)

        for arrow_coord in arrows_cordinates:
            PyMOL_Widget.create_directed_arrows(atom1=arrow_coord[0], atom2=arrow_coord[1],
                                                radius=0.05, name='pairNet',
                                                gap=0.4, hradius=0.4, hlength=0.8,
                                                color='orange')

        for arrow_coord in shortest_path_arrow_coords:
            PyMOL_Widget.create_directed_arrows(atom1=arrow_coord[0], atom2=arrow_coord[1],
                                                radius=0.055,
                                                gap=0.4, hradius=0.4, hlength=0.8,
                                                color='magenta', shortest_path=True)

        for node in processed_path:
            pymol_selection = PyMOL_Widget._build_residue_selection_query(node)
            if pymol_selection is None:
                continue
            PyMOL_Widget.resi_label_add(pymol_selection)

        # MAKE PYMOL VISUALIZATION BETTER
        PyMOL_Widget._pymol.cmd.set('cartoon_oval_length', 0.8)  # default is 1.20)
        PyMOL_Widget._pymol.cmd.set('cartoon_oval_width', 0.2)
        PyMOL_Widget._pymol.cmd.center(selection="all", state=0, origin=1, animate=0)
        PyMOL_Widget._pymol.cmd.zoom('all', buffer=0.0, state=0, complete=0)
        PyMOL_Widget.update()
        PyMOL_Widget.show()

    def get_conservation_scores(self):
        try:
            conserv_pdb_id = self.conservation_PDB_ID_lineEdit.text()
            conser_pdb_chain_id = self.conservation_pdb_chain_id_lineedit.text()
            res_IDs, con_scores = get_conservation_scores(pdb_id=conserv_pdb_id,
                                                          chain_id=conser_pdb_chain_id,
                                                          cutoff=self.conserv_score_doubleSpinBox.value(),
                                                          bound_pdb=self.boundForm_pdb_lineedit.text())

            numrows = len(res_IDs)  # 6 rows in your example
            numcols = 2  # 3 columns in your example
            # Set colums and rows in QTableWidget
            self.residues_conservation_tableWidget.setColumnCount(numcols)
            self.residues_conservation_tableWidget.setRowCount(numrows)

            # Loops to add values into QTableWidget
            for row in range(numrows):
                self.residues_conservation_tableWidget.setItem(row, 0, QTableWidgetItem((res_IDs[row])))
                self.residues_conservation_tableWidget.setItem(row, 1, QTableWidgetItem((str(con_scores[row]))))
        except Exception as Err:
            print("Conservation Score Listing Problem \n", Err)

    def browse_responseTimeFile(self):
        """
            The function provides Main GUI / Upload button activity for select response time file indicated by the user
        """
        try:
            d = file_dialog(self)
            if self.Output_Folder_textEdit.toPlainText().strip() != '' and os.path.exists(
                    self.Output_Folder_textEdit.toPlainText().strip()):
                d.directory = self.Output_Folder_textEdit.toPlainText().strip()
            else:
                d.directory = os.getcwd()
            d.filters = ['Response_Time_Files (*.csv)', 'Tüm Dosyalar (*.*)']
            d.default_filter_index = 0
            d.exec(load=True)

            if d.filename[2] == 'csv':
                self.response_time_lineEdit.setText(d.path)
                return True, d.path

            else:
                Message_Boxes.Critical_message(self, "Error", "This is not a valid response time file",
                                               Style.MessageBox_stylesheet)

        except Exception as exp:
            Message_Boxes.Warning_message(self, "Fatal Error!", str(exp), Style.MessageBox_stylesheet)

    def browse_bound_form_pdbFile(self):
        """
            The function provides Main GUI / Upload button activity for select pdb file indicated by the user
        """
        try:
            d = file_dialog(self)
            if self.Output_Folder_textEdit.toPlainText().strip() != '' and os.path.exists(
                    self.Output_Folder_textEdit.toPlainText().strip()):
                d.directory = self.Output_Folder_textEdit.toPlainText().strip()
            else:
                d.directory = os.getcwd()

            d.filters = ['pdb Files (*.pdb)', 'Tüm Dosyalar (*.*)']
            d.default_filter_index = 0
            d.exec(load=True)

            if d.filename[2] == 'pdb':
                self.boundForm_pdb_lineedit.setText(d.path)
                return True, d.path

            else:
                Message_Boxes.Critical_message(self, "Error", "this is not a pdb file", Style.MessageBox_stylesheet)

        except Exception as exp:
            Message_Boxes.Warning_message(self, "Fatal Error!", str(exp), Style.MessageBox_stylesheet)

    def analysis_output_directory(self):
        try:
            options = QFileDialog.Options()
            options |= QFileDialog.DontUseNativeDialog
            output_file = QFileDialog.getExistingDirectory(options=options)
            self.net_output_directory_lineedit.setText(output_file)
            return True

        except Exception as ins:
            return False

    def node_threshold_use(self):
        if self.node_threshold_checkBox.isChecked():
            self.node_threshold_spinBox.setEnabled(False)
        else:
            self.node_threshold_spinBox.setEnabled(True)

    @staticmethod
    def All_Residues_as_target_Changed(self):
        if not self.all_targets_checkBox.isChecked():
            self.target_res_comboBox.setEnabled(True)
            self.residues_conservation_tableWidget.setEnabled(True)

        if self.all_targets_checkBox.isChecked():
            self.target_res_comboBox.setEnabled(False)
            self.residues_conservation_tableWidget.setEnabled(False)

    @staticmethod
    def All_CPU_Usage_State(self):
        if not self.All_CPU_checkBox.isChecked():
            self.Number_CPU_spinBox_2.setEnabled(True)

        if self.All_CPU_checkBox.isChecked():
            self.Number_CPU_spinBox_2.setEnabled(False)

    # ########################################### ANALYSIS WINDOW FUNCTIONS ############################################

    # ######################################### PERTURBATION WINDOW FUNCTIONS ##########################################
    def maximum_thread_of_system(self):
        # self.Number_CPU_spinBox.setMaximum(mp.cpu_count())
        self.Number_of_thread_for_network_spinBox.setMaximum(mp.cpu_count())

    def precision_combobox_settings(self, eq_md_indexes, per_md_indexes=None):

        if eq_md_indexes is not None:
            if str(self.equ_platform_comboBox.currentText()) == 'CPU':
                self.eq_precision_comboBox.model().item(eq_md_indexes['single']).setEnabled(False)
                self.eq_precision_comboBox.model().item(eq_md_indexes['mixed']).setEnabled(True)
                self.eq_precision_comboBox.model().item(eq_md_indexes['double']).setEnabled(False)
                self.eq_precision_comboBox.setCurrentIndex(eq_md_indexes['mixed'])
                self.All_CPU_checkBox.setEnabled(True)
                self.Number_CPU_spinBox_2.setEnabled(True)
                # self.label_6.setEnabled(True)

            if str(self.equ_platform_comboBox.currentText()) == 'Reference':
                self.eq_precision_comboBox.model().item(eq_md_indexes['single']).setEnabled(False)
                self.eq_precision_comboBox.model().item(eq_md_indexes['mixed']).setEnabled(False)
                self.eq_precision_comboBox.model().item(eq_md_indexes['double']).setEnabled(True)
                self.eq_precision_comboBox.setCurrentIndex(eq_md_indexes['double'])
                self.All_CPU_checkBox.setEnabled(False)
                self.Number_CPU_spinBox_2.setEnabled(False)
                # self.label_6.setEnabled(False)

            if str(self.equ_platform_comboBox.currentText()) in ['CUDA', 'OpenCL']:
                self.eq_precision_comboBox.model().item(eq_md_indexes['single']).setEnabled(True)
                self.eq_precision_comboBox.model().item(eq_md_indexes['mixed']).setEnabled(True)
                self.eq_precision_comboBox.model().item(eq_md_indexes['double']).setEnabled(True)
                self.eq_precision_comboBox.setCurrentIndex(eq_md_indexes['single'])
                self.All_CPU_checkBox.setEnabled(False)
                self.Number_CPU_spinBox_2.setEnabled(False)

    def output_file(self):
        try:
            options = QFileDialog.Options()
            options |= QFileDialog.DontUseNativeDialog
            output_file = QFileDialog.getExistingDirectory(options=options)
            self.Output_Folder_textEdit.setText(output_file)
            return True

        except Exception as ins:
            return False

    def browse_pdbFile(self):
        """
            The function provides Main GUI / Upload button activity for select pdb file indicated by the user
        """
        try:
            fileFormats = "*.pdb *.gro *.psf"  # without commas
            self.pdb_filename, _ = QtWidgets.QFileDialog.getOpenFileName(parent=None,
                                                                         caption="Set pdb file",
                                                                         filter=fileFormats)
            # options = QFileDialog.Options()
            # options |= QFileDialog.DontUseNativeDialog
            # self.pdb_filename, _ = QFileDialog.getOpenFileName(self, "Show The *pdb File", str(os.getcwd()),
            #                                                    "pdb Files (*.pdb)", str(options))

            self.pdb_path = self.pdb_filename
            self.pdb_filename = os.path.splitext(os.path.basename(self.pdb_filename))

            if self.pdb_filename[1] == '.pdb':
                return True, self.pdb_path

            elif self.pdb_filename[1] != "":
                Message_Boxes.Critical_message(self, "Error", "this is not a pdb file", Style.MessageBox_stylesheet)

        except Exception as exp:
            Message_Boxes.Warning_message(self, "Fatal Error!", str(exp), Style.MessageBox_stylesheet)

    def load_sample_for_simulation(self):
        try:
            from pathlib import Path
            current_path = os.path.dirname(os.path.realpath(__file__))
            example_sim_path = os.path.join(
                Path(current_path).parent,
                'example',
                'simulation_demo',
                '2j0w_example.pdb',
            )
            download_sim_path = os.path.join(Path(current_path).parent, 'Download', '2j0w_example.pdb')
            path = example_sim_path if os.path.exists(example_sim_path) else download_sim_path
            output_directory = os.path.join(Path(current_path).parent, 'output')

            if os.path.exists(path):
                self.upload_pdb_lineEdit.setText(path)
                self._pending_res1_index = 342
                self.upload_pdb_from_local(manuel=False)
                try:
                    os.mkdir(output_directory)
                    print("Directory ", output_directory, " Created ")
                except FileExistsError:
                    print("Directory ", output_directory, " already exists")
                if not os.path.exists(output_directory):
                    os.mkdir(output_directory)
                    print("Directory ", output_directory, " Created ")
                else:
                    print("Directory ", output_directory, " already exists")
                self.Output_Folder_textEdit.setPlainText("")

                self.PDB_ID_lineEdit.setText("")
                self.Output_Folder_textEdit.setPlainText(output_directory)

                if self.selected_residues_listWidget.count() == 0:
                    self.selected_residues_listWidget.addItem('SER345A')
                    self.res1_comboBox.setCurrentIndex(342)
                else:
                    self.selected_residues_listWidget.clear()
                    self.selected_residues_listWidget.addItem('SER345A')
                    self.target_res_comboBox.setCurrentIndex(342)

                self.R_factor_ComboBox.setCurrentText("4")
                self.run_duration_spinBox.blockSignals(True)
                self.run_duration_doubleSpinBox.blockSignals(True)
                self.Number_of_steps_spinBox.blockSignals(True)
                self.long_simulation_time_unit.blockSignals(True)
                # self.run_duration_spinBox.setValue(int(1000))
                # self.run_duration_doubleSpinBox.setValue(float(600.000))
                # self.Number_of_steps_spinBox.setValue(int(300000))
                self.run_duration_spinBox.setValue(int(1000))
                self.run_duration_doubleSpinBox.setValue(float(1.000))
                self.Number_of_steps_spinBox.setValue(int(500))

                self.long_simulation_time_unit.setCurrentIndex(1)
                self.run_duration_spinBox.blockSignals(False)
                self.run_duration_doubleSpinBox.blockSignals(False)
                self.Number_of_steps_spinBox.blockSignals(False)
                self.long_simulation_time_unit.blockSignals(False)

            elif not os.path.exists(path):
                answer = Message_Boxes.Information_message(self, "Test PDB File Not Found",
                                                           "The system couldn't find the test PDB file at the "
                                                           "specified location.", Style.MessageBox_stylesheet)

                if answer == QMessageBox.Ok:
                    print("Ok")

        except Exception as Error_info:
            print(Error_info)

    def load_sample_for_analysis(self):
        from pathlib import Path

        def _has_recalculation_inputs(response_csv_path):
            response_dir = os.path.dirname(response_csv_path)
            response_name = os.path.basename(response_csv_path)
            suffix = ""
            if response_name.startswith('responseTimes_') and response_name.lower().endswith('.csv'):
                suffix = response_name[len('responseTimes_'):-4]

            reference_candidates = [os.path.join(response_dir, 'reference_energy_file.csv')]
            modified_candidates = [os.path.join(response_dir, 'modified_energy_file.csv')]
            if suffix:
                reference_candidates.insert(0, os.path.join(response_dir, f'reference_energy_file_{suffix}.csv'))
                modified_candidates.insert(0, os.path.join(response_dir, f'modified_energy_file_{suffix}.csv'))

            has_reference = any(os.path.isfile(candidate) for candidate in reference_candidates)
            has_modified = any(os.path.isfile(candidate) for candidate in modified_candidates)
            return has_reference and has_modified

        current_path = os.path.dirname(os.path.realpath(__file__))
        project_root = Path(current_path).parent.parent
        example_analysis_dir = os.path.join(Path(current_path).parent, 'example', 'analysis_demo')
        preferred_output_dir = os.path.join(
            Path(current_path).parent,
            'output',
            '2j0w_example_fixed_ph7.4_2026-03-15_00-08-18',
        )
        sample_pairs = [
            (
                os.path.join(example_analysis_dir, 'last.pdb'),
                os.path.join(example_analysis_dir, 'responseTimes_4.csv'),
            ),
            (
                os.path.join(example_analysis_dir, 'minimized.pdb'),
                os.path.join(example_analysis_dir, 'responseTimes_4.csv'),
            ),
            (
                os.path.join(preferred_output_dir, 'last.pdb'),
                os.path.join(preferred_output_dir, 'responseTimes_4.csv'),
            ),
            (
                os.path.join(preferred_output_dir, 'minimized.pdb'),
                os.path.join(preferred_output_dir, 'responseTimes_4.csv'),
            ),
            (
                os.path.join(Path(current_path).parent, 'Download', '2j0w_example_fixed_ph7.4.pdb'),
                os.path.join(preferred_output_dir, 'responseTimes_4.csv'),
            ),
            (
                os.path.join(project_root, 'tests', '2DH3', '2DH3.pdb'),
                os.path.join(project_root, 'tests', '2DH3', 'responseTimes_4.csv'),
            ),
        ]

        topology_path = None
        reTime_File = None
        fallback_topology_path = None
        fallback_response_file = None
        for candidate_topology, candidate_response in sample_pairs:
            if os.path.exists(candidate_topology) and os.path.exists(candidate_response):
                if fallback_topology_path is None:
                    fallback_topology_path = candidate_topology
                    fallback_response_file = candidate_response

                if _has_recalculation_inputs(candidate_response):
                    topology_path = candidate_topology
                    reTime_File = candidate_response
                    break

        if topology_path is None and fallback_topology_path is not None:
            topology_path = fallback_topology_path
            reTime_File = fallback_response_file

        if topology_path is None or reTime_File is None:
            Message_Boxes.Warning_message(
                self,
                "Analysis sample files are missing",
                "No valid analysis sample topology/response-time files were found in example, output, Download or tests folders.",
                Style.MessageBox_stylesheet,
            )
            return

        output_directory = os.path.join(Path(current_path).parent, 'output')

        if os.path.exists(topology_path):
            self.boundForm_pdb_lineedit.setText(topology_path)
            self.response_time_lineEdit.setText(reTime_File)
            self._pending_target_index = 255
            self._pending_source_index = 342
            self.upload_boundForm_pdb_from_local(manuel=False)
            try:
                os.mkdir(output_directory)
                print("Directory ", output_directory, " Created ")
            except FileExistsError:
                print("Directory ", output_directory, " already exists")
            if not os.path.exists(output_directory):
                os.mkdir(output_directory)
                print("Directory ", output_directory, " Created ")
            else:
                print("Directory ", output_directory, " already exists")
            self.net_output_directory_lineedit.setText("")
            self.net_output_directory_lineedit.setText(output_directory)

            if self.selected_target_residues_listWidget.count() == 0:
                self.selected_target_residues_listWidget.addItems(
                    ['THR221A', 'PRO231A', 'LYS257A', 'VAL258A'])
                self.target_res_comboBox.setCurrentIndex(255)
            else:
                self.selected_target_residues_listWidget.clear()
                self.selected_target_residues_listWidget.addItems(
                    ['THR221A', 'PRO231A', 'LYS257A', 'VAL258A'])
                self.target_res_comboBox.setCurrentIndex(255)

            self.source_res_comboBox.setCurrentIndex(342)

            self.network_cutoff_spinBox.setValue(float(7.00))
            self.node_threshold_checkBox.setChecked(True)
            self.conservation_PDB_ID_lineEdit.setText("")
            self.conservation_pdb_chain_id_lineedit.setText("")
            self.conserv_score_doubleSpinBox.setValue(float(1.0))
            self.show_analysis_window()

            if not _has_recalculation_inputs(reTime_File):
                Message_Boxes.Information_message(
                    self,
                    "Sample loaded",
                    "Analysis sample loaded successfully.\n"
                    "Recalculate Response Time requires reference_energy_file*.csv and modified_energy_file*.csv, "
                    "which are not available for this sample.",
                    Style.MessageBox_stylesheet,
                )

            # self.run_duration_spinBox.blockSignals(True)
            # self.run_duration_doubleSpinBox.blockSignals(True)
            # self.Number_of_steps_spinBox.blockSignals(True)
            # self.run_duration_doubleSpinBox.setValue(float(1.000))
            # self.Number_of_steps_spinBox.setValue(int(500000))
            # self.run_duration_spinBox.blockSignals(False)
            # self.run_duration_doubleSpinBox.blockSignals(False)
            # self.Number_of_steps_spinBox.blockSignals(False)

    def export_workspace(self):
        try:
            # default_dir = os.path.join(os.path.join(os.path.expanduser('~')), 'Desktop')
            options = QFileDialog.Options()
            options |= QFileDialog.DontUseNativeDialog
            output_file_directory, _ = QFileDialog.getSaveFileName(parent=None, caption="Save configuration file",
                                                                   filter="Config File Name (*.yaml)", options=options)
            template = config_template()

            # INPUTS
            template['Inputs']['topology'] = self.upload_pdb_lineEdit.text().strip()
            template['Inputs']['output directory'] = self.Output_Folder_textEdit.toPlainText().strip()
            template['Inputs']['run duration'] = float(self.run_duration_doubleSpinBox.value())
            template['Inputs']['run duration unit'] = self.long_simulation_time_unit.currentText()
            template['Inputs']['threading number'] = int(self.Number_CPU_spinBox_2.value())
            template['Inputs']['all cpu is active'] = self.All_CPU_checkBox.isChecked()
            template['Inputs']['r factor'] = self.R_factor_ComboBox.currentText().strip()
            template['Inputs']['selected residue'] = [self.selected_residues_listWidget.item(x).text() for x in
                                                      range(self.selected_residues_listWidget.count())]
            template['Inputs']['perturbation run duration'] = int(self.run_duration_spinBox.value())
            template['Inputs']['perturbation run duration unit'] = self.perturb_time_unit_comboBox.currentText()
            template['Inputs']['perturbation threading number'] = int(self.Number_CPU_spinBox.value())
            template['Inputs']['perturbation all cpu is active'] = self.perturbation_All_CPU_checkBox.isChecked()

            # SIMULATION
            template['Simulation']['equilibrium platform'] = self.equ_platform_comboBox.currentText()
            template['Simulation']['equilibrium precision'] = self.eq_precision_comboBox.currentText()
            template['Simulation']['equilibrium device id'] = self.Device_Number_comboBox.currentText()
            template['Simulation']['equilibrium use device id'] = self.Device_ID_checkBox.isChecked()
            template['Simulation']['perturbation platform'] = self.per_platform_comboBox.currentText()
            template['Simulation']['perturbation precision'] = self.per_precision_comboBox.currentText()
            template['Simulation']['perturbation device id'] = self.per_Device_Number_lineEdit.text().strip()
            template['Simulation']['perturbation use device id'] = self.per_Device_ID_checkBox.isChecked()
            template['Simulation']['protein forcefield'] = self.protein_forcefield_comboBox.currentText()
            template['Simulation']['water forcefield'] = self.water_forcefield_comboBox.currentText()
            template['Simulation']['water geometry padding'] = self.water_padding_lineEdit.text().strip()
            template['Simulation']['equilibrium integrator'] = self.integrator_kind_comboBox.currentText()
            template['Simulation']['equilibrium time step'] = self.integrator_time_step_lineEdit.text().strip()
            template['Simulation']['equilibrium time step unit'] = self.integrator_time_step_unit.currentText()
            template['Simulation'][
                'equilibrium additional integrator parameters'] = self.Additional_Integrators_checkBox.isChecked()
            template['Simulation']['friction'] = self.friction_lineEdit.text().strip()
            template['Simulation']['temperature'] = self.temperature_lineEdit.text().strip()
            template['Simulation']['nonbonded method'] = self.nonBounded_Method_comboBox.currentText()
            template['Simulation']['constraints'] = self.system_constraints_comboBox.currentText()
            template['Simulation']['rigid water is active'] = self.rigid_water_checkBox.isChecked()
            template['Simulation']['nonbonded cutoff'] = self.nonbounded_CutOff_lineEdit.text()
            template['Simulation']['switching distance'] = self.switching_distance_lineEdit.text()
            template['Simulation']['use switching distance'] = self.use_switching_checkBox.isChecked()
            template['Simulation']['number of simulation steps'] = self.Number_of_steps_spinBox.value()
            template['Simulation']['minimize'] = self.minimize_checkBox.isChecked()
            template['Simulation']['minimization max iterations'] = self.Max_minimize_iter_lineEdit.text()
            template['Simulation']['Equilibrate'] = self.equilubrate_checkBox.isChecked()
            template['Simulation']['Equilibration steps'] = self.Max_equilubrate_steps_textEdit.toPlainText()
            template['Simulation']['DCD report'] = self.DCD_Reporter_checkBox.isChecked()
            template['Simulation']['XTC report'] = self.XTC_Reporter_checkBox.isChecked()
            template['Simulation']['state data report'] = self.State_Data_Reporter_checkBox.isChecked()
            template['Simulation']['DCD writing frequency'] = self.DCD_write_freq_lineEdit.text()
            template['Simulation']['DCD output name'] = self.DCD_Output_Name_lineEdit.text()
            template['Simulation']['XTC writing frequency'] = self.XTC_write_freq_lineEdit.text()
            template['Simulation']['XTC output name'] = self.XTC_Output_Name_lineEdit.text()
            template['Simulation']['state data frequency'] = self.StateData_frequency_lineEdit.text()

            write_output_configuration_file(file_path=output_file_directory, template_yml=template)

        except Exception as err:
            pass

    def import_workspace(self):
        try:
            config_file_path, _ = QtWidgets.QFileDialog.getOpenFileName(parent=None,
                                                                        caption="Select configuration file to load",
                                                                        filter="*.yaml")
            template = read_output_configuration_file(file_path=config_file_path)

            if os.path.exists(template['Inputs']['topology']):
                self.upload_pdb_lineEdit.setText(template['Inputs']['topology'])
                self.upload_pdb_from_local(manuel=False)

            try:
                os.mkdir(template['Inputs']['output directory'])
                self.Output_Folder_textEdit.setPlainText(template['Inputs']['output directory'])
            except FileExistsError:
                self.Output_Folder_textEdit.setPlainText(template['Inputs']['output directory'])

            self.run_duration_doubleSpinBox.setValue(float(template['Inputs']['run duration']))
            self.long_simulation_time_unit.setCurrentIndex(
                self.long_simulation_time_unit.findText(
                    str(template['Inputs']['run duration unit']), QtCore.Qt.MatchFixedString))
            self.Number_CPU_spinBox_2.setValue(int(template['Inputs']['threading number']))
            self.All_CPU_checkBox.setChecked(template['Inputs']['all cpu is active'])
            self.R_factor_ComboBox.setCurrentText(str(template['Inputs']['r factor']))
            self.selected_residues_listWidget.addItems(template['Inputs']['selected residue'])
            self.run_duration_spinBox.setValue(int(template['Inputs']['perturbation run duration']))
            self.perturb_time_unit_comboBox.setCurrentIndex(
                self.perturb_time_unit_comboBox.findText(
                    str(template['Inputs']['perturbation run duration unit']), QtCore.Qt.MatchFixedString))
            self.Number_CPU_spinBox.setValue(int(template['Inputs']['perturbation threading number']))
            self.perturbation_All_CPU_checkBox.setChecked(template['Inputs']['perturbation all cpu is active'])

        except Exception as Err:
            pass

    @staticmethod
    def PDB_ID_lineEdit(self):
        """
             The function provides enable or disable of Upload Button according to the entry of user
        """
        self.text_edit = self.PDB_ID_lineEdit.text()
        if self.text_edit == "":
            self.upload_pdb_Button.setEnabled(True)
            self.upload_pdb_lineEdit.setEnabled(True)
            self.label_32.setEnabled(True)
        else:
            self.upload_pdb_Button.setEnabled(False)
            self.upload_pdb_lineEdit.setEnabled(False)
            self.label_32.setEnabled(False)

    @staticmethod
    def Fetch_PDB_File(self):
        global selected_chains, fetched_pdb
        try:
            fetch_pdb_ID = self.PDB_ID_lineEdit.text()
            if len(fetch_pdb_ID) == 4:

                fetch_result = pdb_Tools.fetch_pdb(self, fetch_pdb_ID)
                if fetch_result != False:

                    if path.exists(fetch_result):
                        fixer = PDBFixer(fetch_result)
                        fixer.removeHeterogens(keepWater=False)

                        modeller = Modeller(fixer.topology, fixer.positions)
                        chains = [r.id for r in modeller.topology.chains()]

                        checked_list = ChecklistDialog('Select the chain (s) to be used in the system', chains,
                                                       checked=True)

                        pdb_fix_dialog_answer = checked_list.exec_()

                        if pdb_fix_dialog_answer == QDialog.Accepted:
                            selected_chains = [str(s) for s in checked_list.choices]
                            delete_chains = list(set(chains) - set(selected_chains))
                            fetched_pdb = pdb_Tools.fetched_pdb_fix(self, fetch_result,
                                                                    self.Output_Folder_textEdit.toPlainText(),
                                                                    ph=self.pH_doubleSpinBox.value(),
                                                                    chains_to_remove=delete_chains)

                            self.upload_pdb_lineEdit.setText(fetched_pdb)
                            self.combobox = Helper_Functions.fill_residue_combobox(self, fetched_pdb)
                            for i in self.combobox:
                                self.res1_comboBox.addItem(str(i))
                            self.res1_comboBox.clear()  # delete all items from comboBox
                            self.res1_comboBox.addItems(self.combobox)  # add the actual content of self.comboData
                            InputFile(fetch_result)
                            return fetched_pdb

                        elif pdb_fix_dialog_answer == QDialog.Rejected:
                            modified_pdb = pdb_Tools.fetched_pdb_fix(self, fetch_result,
                                                                     self.Output_Folder_textEdit.toPlainText(),
                                                                     ph=self.pH_doubleSpinBox.value(),
                                                                     chains_to_remove=None)

                            self.upload_pdb_lineEdit.setText(modified_pdb)

                            self.combobox = Helper_Functions.fill_residue_combobox(self, modified_pdb)
                            for i in self.combobox:
                                self.res1_comboBox.addItem(str(i))
                            self.res1_comboBox.clear()  # delete all items from comboBox
                            self.res1_comboBox.addItems(self.combobox)  # add the actual content of self.comboData

                            InputFile(modified_pdb)
                            return modified_pdb

                        return None

                else:
                    return None

            if len(fetch_pdb_ID) != 4:
                Message_Boxes.Information_message(self, 'Wrong pdb id', 'PDB ID should be provided as 4 letters',
                                                  Style.MessageBox_stylesheet)
                return False

        except Exception as instance:
            Message_Boxes.Critical_message(self, 'An error occurred while fetching the pdb file.', repr(instance),
                                           Style.MessageBox_stylesheet)

    @staticmethod
    def Equilibration_On_Off_Changed(self):
        if not self.equilubrate_checkBox.isChecked():
            self.MD_Run_Time_GroupBox.setEnabled(False)
            self.groupBox_18.setEnabled(False)

        if self.equilubrate_checkBox.isChecked():
            self.MD_Run_Time_GroupBox.setEnabled(True)
            self.groupBox_18.setEnabled(True)

    @staticmethod
    def UseDeviceID_On_Off_Changed(self):
        if not self.Device_ID_checkBox.isChecked():
            self.Device_Number_comboBox.setEnabled(False)

        if self.Device_ID_checkBox.isChecked():
            self.Device_Number_comboBox.setEnabled(True)

    @staticmethod
    def minimize_Step_isVisible(self):
        if not self.minimize_checkBox.isChecked():
            self.Max_minimize_iter_lineEdit.setEnabled(False)
        else:
            self.Max_minimize_iter_lineEdit.setEnabled(True)

    @staticmethod
    def DCD_Reporter_Changed(self):
        if not self.DCD_Reporter_checkBox.isChecked():
            self.DCD_Reporter_Options_groupBox.setEnabled(False)

        if self.DCD_Reporter_checkBox.isChecked():
            self.DCD_Reporter_Options_groupBox.setEnabled(True)

    @staticmethod
    def XTC_Reporter_Changed(self):
        if not self.XTC_Reporter_checkBox.isChecked():
            self.XTC_Reporter_Options_groupBox.setEnabled(False)

        if self.XTC_Reporter_checkBox.isChecked():
            self.XTC_Reporter_Options_groupBox.setEnabled(True)

    @staticmethod
    def State_Data_Reporter_Changed(self):
        if not self.State_Data_Reporter_checkBox.isChecked():
            self.State_Data_Reporter_Options_groupBox.setEnabled(False)

        if self.State_Data_Reporter_checkBox.isChecked():
            self.State_Data_Reporter_Options_groupBox.setEnabled(True)

    @staticmethod
    def Send_Available_Platforms_to_GUI(self):
        self.platforms, self.plt_speeds = Helper_Functions.available_platforms(self)
        self.platform_list_on_the_program = [self.equ_platform_comboBox.itemText(i) for i in
                                             range(self.equ_platform_comboBox.count())]

        for item_no, i in enumerate(self.platform_list_on_the_program):
            if i not in self.platforms:
                self.equ_platform_comboBox.model().item(int(item_no)).setEnabled(False)
                self.equ_platform_comboBox.setCurrentIndex(item_no + 1)

            if i in self.platforms:
                self.equ_platform_comboBox.setItemData(item_no, str("Estimated Speed For This Devices Is "
                                                                    + str(self.plt_speeds[i])), Qt.ToolTipRole)

    @staticmethod
    def platform_comboBox_Changed(self, eq_md, per_md):
        if eq_md:
            if self.equ_platform_comboBox.currentText() in ["CPU", "Reference"]:
                self.Device_Number_comboBox.setEnabled(False)
                self.Device_ID_checkBox.setEnabled(False)

            else:
                self.Device_Number_comboBox.setEnabled(True)
                self.Device_ID_checkBox.setEnabled(True)

        if per_md:
            if self.per_platform_comboBox.currentText() in ["CPU", "Reference"]:
                self.per_Device_Number_lineEdit.setEnabled(False)
                self.per_Device_ID_checkBox.setEnabled(False)

            else:
                self.per_Device_Number_lineEdit.setEnabled(True)
                self.per_Device_ID_checkBox.setEnabled(True)

    def add_residue_toList(self):
        if str(self.res1_comboBox.currentText()) != "":
            items = []
            for x in range(self.selected_residues_listWidget.count()):
                items.append(self.selected_residues_listWidget.item(x).text())
            if str(self.res1_comboBox.currentText()) not in items:
                self.selected_residues_listWidget.addItem(str(self.res1_comboBox.currentText()))

    def discard_residue_fromList(self):
        listItems = self.selected_residues_listWidget.selectedItems()
        if not listItems:
            return
        for item in listItems:
            self.selected_residues_listWidget.takeItem(self.selected_residues_listWidget.row(item))

    def add_residue_to_target_List(self):
        if str(self.target_res_comboBox.currentText()) != "":
            items = []
            for x in range(self.selected_target_residues_listWidget.count()):
                items.append(self.selected_target_residues_listWidget.item(x).text())
            if str(self.target_res_comboBox.currentText()) not in items:
                self.selected_target_residues_listWidget.addItem(str(self.target_res_comboBox.currentText()))

    def discard_residue_from_target_List(self):
        listItems = self.selected_target_residues_listWidget.selectedItems()
        if not listItems:
            return
        for item in listItems:
            self.selected_target_residues_listWidget.takeItem(self.selected_target_residues_listWidget.row(item))

    def number_of_steps_changed_from_quick(self):
        global new_step

        current_step = self.run_duration_doubleSpinBox.value()
        current_time_unit = self.long_simulation_time_unit.currentText()
        current_integrator_time_step_value = float(self.integrator_time_step_lineEdit.text())

        if current_time_unit == 'nanosecond':
            new_step = int((current_step / current_integrator_time_step_value) * 1000000)

        if current_time_unit == 'picosecond':
            new_step = int((current_step / current_integrator_time_step_value) * 1000)

        self.Number_of_steps_spinBox.setValue(new_step)

    def number_of_steps_changed_from_advanced(self):
        global new_time
        current_step = int(self.Number_of_steps_spinBox.value())  # 1 ns
        current_time_unit = self.long_simulation_time_unit.currentText()  # ns
        current_integrator_time_step_value = float(self.integrator_time_step_lineEdit.text())  # 2 fs

        if current_time_unit == 'nanosecond':
            new_time = float((current_step * current_integrator_time_step_value) / 1000000)

        if current_time_unit == 'picosecond':
            new_time = float((current_step * current_integrator_time_step_value) / 1000)

        self.run_duration_doubleSpinBox.setValue(new_time)

    @staticmethod
    def initialize_advanced_platform_options(main_window):
        """Initialize advanced platform options (Determinism, UseCpuPme, etc.) in UI."""
        try:
            from ._advanced_platform_options import initialize_advanced_platform_options as init_advanced_opts
            init_advanced_opts(main_window)
        except Exception as e:
            print(f"Warning: Could not initialize advanced platform options: {e}")


class InputFile:
    fetch_result = False

    def __init__(self, input_file_geting):
        self.input_file_geting = input_file_geting
        self.pdb_file, fetch_result = self.upload_fetched_pdb(self.input_file_geting)

    def upload_fetched_pdb(self, file_geting):
        if file_geting != None:
            fetch_result = True
            self.pdb_file = file_geting
            return self.pdb_file, fetch_result


def download_pdb_file(pdb_id, compressed, dest_folder):
    # Ensure download folder exists
    try:
        os.makedirs(dest_folder)
    except OSError as e:
        print("Ignore OSError raised if it already exists")

    filename = '%s.pdb' % pdb_id
    # Add .gz extension if compressed
    if compressed:
        filename = '%s.gz' % filename
    url = 'https://files.rcsb.org/download/%s' % filename
    destination_file = os.path.join(dest_folder, filename)
    # Download the file
    urlretrieve(url, destination_file)

    return destination_file


class pdb_Tools:

    def fetch_pdb(self, pdb_id):

        """
        :param pdb_id: 4 letters PDB id provides protein structure file from www.rcsb.org
        :return: unziped pdb file for load
        """
        try:
            from pathlib import Path
            import os
            from prody.proteins.localpdb import fetchPDB, pathPDBFolder
            from multiprocessing import Process
            Download_folder = os.path.join(os.getcwd(), 'Download')

            try:
                os.makedirs(Download_folder)
            except FileExistsError:
                print("directory already exists")
                pass

            fetched_pdb_file = download_pdb_file(pdb_id, compressed=False, dest_folder=Download_folder)
            return fetched_pdb_file

        except Exception as instance:
            Message_Boxes.Critical_message(self, 'An error occurred while fetching the pdb file.', repr(instance),
                                           Style.MessageBox_stylesheet)
            return False
        

    def fetched_pdb_fix(self, file_pathway, output_path=None, ph=7, chains_to_remove=None):
        """
        Args:
            :param file_pathway: Pathway for manipulating your fetched PDB files
            :param chains_to_remove: Selected chains will be deleted
            :param ph: Selected pH value will be applied to the structure's Hydrogens
        Returns:
            :param output_path: The manipulated PDB file will return as full path if specified
                                otherwise will return the already existing path
        """

        def remove_hetatoms_from_pdb(file_content):
            """
            Remove all HETATOM records from the PDB file content.

            Args:
                :param file_content: Content of the input PDB file as a string
            Returns:
                :return: Content of the PDB file with HETATOMs removed as a string
            """
            lines = file_content.splitlines()
            filtered_lines = [line for line in lines if not line.startswith("HETATM")]
            return '\n'.join(filtered_lines)

        # Read the original PDB file
        with open(file_pathway, 'r') as infile:
            pdb_content = infile.read()

        # Create a temporary file for the original PDB content
        with tempfile.NamedTemporaryFile(delete=False, mode='w', suffix='.pdb') as temp_file:
            temp_file.write(pdb_content)
            temp_file_path = temp_file.name

        print("Creating PDBFixer...")
        fixer = PDBFixer(temp_file_path)
        print("Finding missing residues...")

        if chains_to_remove:
            print("toDelete: %s" % chains_to_remove)
            fixer.removeChains(chainIds=chains_to_remove)

        fixer.findMissingResidues()

        chains = list(fixer.topology.chains())
        keys = fixer.missingResidues.keys()
        for key in list(keys):
            chain = chains[key[0]]
            if key[1] == 0 or key[1] == len(list(chain.residues())):
                del fixer.missingResidues[key]

        print("Finding nonstandard residues...")
        fixer.findNonstandardResidues()
        print("Replacing nonstandard residues...")
        fixer.replaceNonstandardResidues()
        print("Removing heterogens...")
        fixer.removeHeterogens(keepWater=False)

        # Remove HETATOM records directly from the PDB content
        cleaned_pdb_content = remove_hetatoms_from_pdb(pdb_content)

        # Create a temporary file for the cleaned PDB content
        with tempfile.NamedTemporaryFile(delete=False, mode='w', suffix='.pdb') as temp_file:
            temp_file.write(cleaned_pdb_content)
            cleaned_pdb_path = temp_file.name

        # Read the cleaned PDB content with HETATOM removed
        fixer = PDBFixer(cleaned_pdb_path)

        if chains_to_remove:
            print("toDelete (cleaned): %s" % chains_to_remove)
            fixer.removeChains(chainIds=chains_to_remove)

        print("Adding missing atoms...")
        fixer.findMissingResidues()
        print("Adding missing hydrogens...")
        fixer.addMissingHydrogens(pH=ph)

        print("Writing PDB file...")

        if output_path:
            path = os.path.join(output_path, f"{os.path.basename(file_pathway).split('.')[0]}_fixed_ph{ph}.pdb")
            with open(path, "w") as outfile:
                PDBFile.writeFile(fixer.topology, fixer.positions, outfile, keepIds=True)
            return path

        if not output_path:
            new_outpath = os.path.join(os.path.dirname(file_pathway),
                                       f"{os.path.basename(file_pathway).split('.')[0]}_fixed_ph{ph}.pdb")
            with open(new_outpath, "w") as outfile:
                PDBFile.writeFile(fixer.topology, fixer.positions, outfile, keepIds=True)
            return new_outpath

    @staticmethod
    def initialize_advanced_platform_options(main_window):
        """Initialize advanced platform options (Determinism, UseCpuPme, etc.) in UI."""
        try:
            from ._advanced_platform_options import initialize_advanced_platform_options as init_advanced_opts
            init_advanced_opts(main_window)
        except Exception as e:
            print(f"Warning: Could not initialize advanced platform options: {e}")