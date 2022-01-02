import os.path

import networkx as nx
from PySide2.QtWidgets import QFileDialog, QDialog, QTableWidgetItem
from PySide2.QtCore import Qt, Slot
import csv
from pathlib import Path
from pdbfixer import PDBFixer
from os import path
from urllib.request import urlretrieve
from openmm.app import *
from .checkBox_menu import *
from .message import Message_Boxes
from .PyMolWidget import PymolQtWidget
import multiprocessing as mp
from analysis.pdbsum_conservation_puller import get_conservation_scores
from analysis.createRNetwork import (Multi_Task_Engine, intersection_of_directed_networks, Pymol_Visualize_Path,
                                     Shortest_Path_Visualize)
from .config import write_output_configuration_file, read_output_configuration_file, config_template
from ui_main import *


class Helper_Functions():

    def fill_residue_combobox(self, pdb_path):
        from prody.proteins.pdbfile import parsePDB
        combobox = []

        pdb = parsePDB(pdb_path)
        protein = pdb.select('protein')
        for model in protein.getHierView():
            for chain in model:
                # print(chain)
                # self.combobox.append(str(model).split(" ")[1] + str(chain))
                combobox.append(str(chain).replace(" ", "") + str(model).split(" ")[1])

                # self.combobox_2.append(str(model).split(" ")[1] + str(chain))
        return combobox

    def available_platforms(self):
        """
            :return: The function will return available platforms on your system for OpenMM Engine
        """
        import openmm
        avail_plt_and_speed = dict()

        for index in range(openmm.Platform.getNumPlatforms()):
            avail_plt_and_speed[
                (openmm.Platform.getPlatform(index).getName())] = openmm.Platform.getPlatform(
                index).getSpeed()

        return avail_plt_and_speed.keys(), avail_plt_and_speed


class Functions(MainWindow):
    global active_workers

    # ########################################### ANALYSIS WINDOW FUNCTIONS ############################################

    @Slot()
    def on_started(self):
        if self.active_workers == 0:
            self.active_workers += 1
            self.progress = QProgressDialog('Work in progress...', None, 0, 0, self)
            self.progress.setWindowTitle("Calculation")
            self.progress.setWindowModality(Qt.WindowModal)
            self.progress.show()
            self.progress.setValue(0)

        else:
            self.active_workers += 1

    @Slot()
    def progress_fn(self, progress_on_net_calc):
        print(progress_on_net_calc, "IS DONE")
        if self.active_workers == 0:
            print("NOW FINISHED ON PROGRESS")

    def thread_complete(self):
        self.active_workers -= 1
        if self.active_workers == 0 and len(self.network_holder) != 0:
            print("===============================")
            print(self.log_holder, self.network_holder)
            print("===============================")
            self.plot_signal.plot_network.emit()
            self.progress.cancel()

    def print_output(self, s):
        self.network_holder.append(s[0])
        self.log_holder.append(s[1])

    def calculate_intersection_network(self):
        global intersection_graph, output_folder_directory, network_holder

        self.active_workers = 0
        self.network_holder = []
        self.log_holder = []

        self.number_of_threads = self.Number_of_thread_for_network_spinBox.value()
        self.pdb = self.boundForm_pdb_lineedit.text()  # PDB file path --> "BOUND FORM OF STRUCTURE"
        self.cutoff = self.network_cutoff_spinBox.value()  # Add edges between nodes if they are within cutoff range
        self.retime_file = self.response_time_lineEdit.text()  # Response time file path
        self.outputFileName = self.PPI_Network_name_lineedit.text()  # Protein general graph according to cut off value
        self.output_directory = self.net_output_directory_lineedit.text()
        self.source = self.source_res_comboBox.currentText()[:-1]  # One of the perturbed residues
        self.node_threshold = self.node_threshold_spinBox.value()  # None or an Integer
        self.node_threshold_use_condition = self.node_threshold_checkBox.isChecked()
        self.all_residue_as_target = self.all_targets_checkBox.isChecked()

        if self.node_threshold_use_condition:
            self.node_threshold = None

        verbose_condition = True  # True or False

        if self.all_residue_as_target:
            target_residues = [self.target_res_comboBox.itemText(i)[:-1] for i in
                               range(self.target_res_comboBox.count())]

        if not self.all_residue_as_target:
            target_residues = [self.selected_target_residues_listWidget.item(x).text()[:-1]
                               for x in range(self.selected_target_residues_listWidget.count())]  # None or residue list

        use_conservation = self.use_conservation_checkBox.isChecked()

        pdb_id = self.conservation_PDB_ID_lineEdit.text()  # FREE OR BOUND FORM OF PDB, like '2EB8'
        chain = self.conservation_pdb_chain_id_lineedit.text()  # FOR PULLING CONS. SCORES INDICATE CHAIN ID OF PDB
        conservation_threshold = self.conserv_score_doubleSpinBox.value()
        save_conservation_scores = False

        self.create_output = True  # Supports True or False Conditions for creation of all networks (*.gml) on a folder

        general_output_folder = os.path.join(self.output_directory, 'network_outputs')
        Path(general_output_folder).mkdir(parents=True, exist_ok=True)

        folder_name = "output_%s" % self.source
        output_folder_directory = os.path.join(general_output_folder, folder_name)
        Path(output_folder_directory).mkdir(parents=True, exist_ok=True)

        engine = Multi_Task_Engine(pdb_file=self.pdb, cutoff=self.cutoff, reTimeFile=self.retime_file,
                                   source=self.source,
                                   node_threshold=self.node_threshold, verbose=verbose_condition,
                                   outputFileName=self.outputFileName, write_outputs=self.create_output,
                                   output_directory=output_folder_directory)

        network, resId_List, len_of_reTimes = engine.calculate_general_network()

        if len(resId_List) == len_of_reTimes:

            if use_conservation:
                res_IDs, con_scores = get_conservation_scores(pdb_id=pdb_id, chain_id=chain,
                                                              cutoff=conservation_threshold, bound_pdb=self.pdb)
                if save_conservation_scores:
                    rows = zip(res_IDs, con_scores)
                    with open(os.path.join(output_folder_directory, 'conservation_%s.csv' % pdb_id), "w",
                              newline='') as f:
                        writer = csv.writer(f)
                        for row in rows:
                            writer.writerow(row)
                intersection_resIDs = set.intersection(set(res_IDs), set(target_residues))

                engine.run_pairNet_calc(intersection_resIDs)
                for work in engine.Work:
                    work.signals.progress_on_net_calc.connect(lambda complete: Functions.progress_fn(self, complete))
                    work.signals.work_started.connect(lambda: Functions.on_started(self))
                    work.signals.result.connect(lambda x: Functions.print_output(self, x))
                    work.signals.finished.connect(lambda: Functions.thread_complete(self))
                    self.threadpool.start(work)

            if not use_conservation:
                engine.run_pairNet_calc(target_residues)

                for work in engine.Work:
                    work.signals.progress_on_net_calc.connect(lambda complete: Functions.progress_fn(self, complete))
                    work.signals.work_started.connect(lambda: Functions.on_started(self))
                    work.signals.result.connect(lambda x: Functions.print_output(self, x))
                    work.signals.finished.connect(lambda: Functions.thread_complete(self))

                    self.threadpool.start(work)

            del engine

        else:
            Message_Boxes.Warning_message(self, "Mismatch Error!", "The number of residues in the topology file you "
                                                                   "have provided is not equal to the response time file.",
                                          Style.MessageBox_stylesheet)
            del engine

    def plot_networks(self):
        print("================***===============")
        print(self.log_holder, self.network_holder)
        print("================***===============")
        clean_graph_list = []

        if self.node_threshold is not None:
            for i in self.network_holder:
                if len(i.nodes()) > self.node_threshold:
                    clean_graph_list.append(i)

        if self.node_threshold is None:
            for i in self.network_holder:
                clean_graph_list.append(i)

        # CREATE AN INTERSECTION GRAPH AND WRITE TO GML FILE
        if len(clean_graph_list) > 0:
            intersection_graph, all_graph_list = intersection_of_directed_networks(clean_graph_list)
            if self.create_output:
                nx.write_gml(intersection_graph, os.path.join(output_folder_directory, 'intersection_graph.gml'))

            # ############################################
            tab = QtWidgets.QWidget()
            tab.setObjectName("Analysis_" + str(self.tab_count_on_analysis))

            self.analysis_TabWidget.tabBar().setTabButton(0, QTabBar.RightSide, None)
            self.tab_count_on_analysis = self.analysis_TabWidget.count()

            horizontalLayout = QtWidgets.QHBoxLayout(tab)
            horizontalLayout.setObjectName("horizontalLayout_" + str(self.tab_count_on_analysis))
            gridLayout = QtWidgets.QGridLayout()
            gridLayout.setObjectName("gridLayout_" + str(self.tab_count_on_analysis))

            label = QtWidgets.QLabel(tab)
            label.setText("All Available Intersection Shortest Path(s)")
            label.setMinimumSize(QtCore.QSize(0, 22))
            label.setMaximumSize(QtCore.QSize(450, 22))
            label.setStyleSheet("QLabel {\n"
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
            label.setObjectName("label_" + str(self.tab_count_on_analysis))
            gridLayout.addWidget(label, 2, 0, 1, 1)

            shortest_path_listWidget = QtWidgets.QListWidget(tab)
            shortest_path_listWidget.setMaximumSize(QtCore.QSize(450, 16777215))
            shortest_path_listWidget.setObjectName("shortest_path_listWidget")
            gridLayout.addWidget(shortest_path_listWidget, 1, 0, 1, 1)

            intersection_path_listWidget = QtWidgets.QListWidget(tab)
            intersection_path_listWidget.setMaximumSize(QtCore.QSize(450, 16777215))
            intersection_path_listWidget.setObjectName("intersection_path_listWidget")
            gridLayout.addWidget(intersection_path_listWidget, 3, 0, 1, 1)

            label_2 = QtWidgets.QLabel(tab)
            label_2.setText("All Available Shortest Path(s)")
            label_2.setMinimumSize(QtCore.QSize(0, 22))
            label_2.setMaximumSize(QtCore.QSize(450, 22))
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
            gridLayout.addWidget(label_2, 0, 0, 1, 1)

            label_3 = QtWidgets.QLabel(tab)
            label_3.setText("Energy Dissipation Curve")
            label_3.setMinimumSize(QtCore.QSize(0, 22))
            label_3.setMaximumSize(QtCore.QSize(450, 22))
            label_3.setStyleSheet("QLabel {\n"
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
            label_3.setObjectName("label_27")
            gridLayout.addWidget(label_3, 4, 0, 1, 1)

            # --> Dissipation Widget
            dissipation_curve_widget = WidgetPlot(self)
            dissipationVerticalLayout = QtWidgets.QVBoxLayout()
            dissipationVerticalLayout.addWidget(dissipation_curve_widget.toolbar)
            dissipationVerticalLayout.addWidget(dissipation_curve_widget.canvas)
            widget = QtWidgets.QWidget()
            widget.setLayout(dissipationVerticalLayout)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(widget.sizePolicy().hasHeightForWidth())
            widget.setSizePolicy(sizePolicy)
            widget.setMinimumSize(QtCore.QSize(0, 300))
            widget.setMaximumSize(QtCore.QSize(450, 350))
            dissipation_curve_widget.setObjectName("dissipation_curve_widget")
            gridLayout.addWidget(widget, 5, 0, 1, 1)

            possible_path = str(self.response_time_lineEdit.text())
            if os.path.exists(possible_path.strip()) and possible_path.split('.')[-1] == 'csv':
                source_residue = self.source_res_comboBox.currentText()
                row, col, Response_Count = getResponseTimeGraph(possible_path)

                if source_residue == '':
                    dissipation_curve_widget.canvas.plot(Response_Count, source_residue=None)
                if source_residue != '':
                    dissipation_curve_widget.canvas.plot(Response_Count, source_residue=source_residue)

            # --> 3D View Frame
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
            gridLayout.addWidget(pyMOL_3D_analysis_frame, 0, 1, 6, 1)
            horizontalLayout.addLayout(gridLayout)
            self.analysis_TabWidget.addTab(tab, "Analysis " + str(self.tab_count_on_analysis))

            # ################################# ==> START - 3D WIDGETS LOCATING <== ################################## #
            Protein3DNetworkView = PymolQtWidget(self)
            verticalLayoutProteinNetworkView = QVBoxLayout(pyMOL_3D_analysis_frame)
            verticalLayoutProteinNetworkView.addWidget(Protein3DNetworkView)
            self.setLayout(verticalLayoutProteinNetworkView)
            Protein3DNetworkView.loadMolFile(self.boundForm_pdb_lineedit.text())
            Protein3DNetworkView.update()
            Protein3DNetworkView.show()
            verticalLayoutProteinNetworkView.setContentsMargins(0, 0, 0, 0)

            # matplotlib_widget = WidgetPlot(self)
            # verticalLayout.addWidget(self.matplotlib_widget.toolbar)
            # verticalLayout.addWidget(self.matplotlib_widget.canvas)

            source_res = self.source_res_comboBox.currentText()[:-1]
            target_res_list = [self.selected_target_residues_listWidget.item(x).text()[:-1]
                               for x in range(self.selected_target_residues_listWidget.count())]

            # #################################################
            shrotest_str_form = ''

            colors = ['#957DAD', '#D291BC', '#FEC8D8', '#8dbdc7', '#B3ABCF', '#b5b1c8', '#e8abb5']

            all_paths = []
            for graph_i in all_graph_list:
                if isinstance(graph_i, nx.classes.digraph.DiGraph):
                    for cnt, target_i in enumerate(target_res_list):
                        try:
                            if nx.has_path(graph_i, source_res, target_i):
                                sp = nx.shortest_path(graph_i, source_res, target_i)
                                for res_id in range(len(sp)):
                                    if res_id == len(sp) - 1:
                                        shrotest_str_form += '%s' % sp[res_id]
                                    else:
                                        shrotest_str_form += '%s --> ' % sp[res_id]

                                item = QListWidgetItem(shrotest_str_form)
                                item.setBackground(QColor(colors[cnt]))
                                shortest_path_listWidget.addItem(item)  # print("SHORTEST PATH: ", sp)
                                shrotest_str_form = ''

                                all_paths.append(
                                    'Source: %s  Target: %s\nTotal node number of source-target pair network is : %s' % (
                                        source_res, target_i, len(graph_i.nodes())))

                        except Exception as err:
                            print("SHORTEST PATH LOG: ", err)

            all_path_string = ''
            for j in all_paths:
                if j != '':
                    all_path_string = all_path_string + '\n' + j
                    print("all: ", all_path_string)
            Message_Boxes.Information_message(self, "DONE !", all_path_string, Style.MessageBox_stylesheet)

            intersect_shrotest_str_form = ''
            for target_i in target_res_list:
                try:
                    if isinstance(graph_i, nx.classes.digraph.DiGraph):
                        isp = nx.shortest_path(intersection_graph, source_res, target_i)
                        for res_id in range(len(isp)):
                            if res_id == len(isp) - 1:
                                intersect_shrotest_str_form += '%s' % isp[res_id]
                            else:
                                intersect_shrotest_str_form += '%s --> ' % isp[res_id]
                        intersection_path_listWidget.addItem(intersect_shrotest_str_form)
                        intersect_shrotest_str_form = ''

                except Exception as err:
                    print("INTERSECTION SHORTEST PATH LOG: ", err)

            shortest_path_listWidget.itemDoubleClicked.connect(
                lambda item: Functions.show_shortest_paths_on_3D_ProteinView(self, item, Protein3DNetworkView))
            intersection_path_listWidget.itemDoubleClicked.connect(
                lambda item: Functions.show_shortest_paths_on_3D_ProteinView(self, item, Protein3DNetworkView))
            # ############################################

            try:
                if isinstance(graph_i, nx.classes.digraph.DiGraph):
                    if len(intersection_graph.nodes()) > 0:
                        try:
                            arrows_cordinates, intersection_node_list = Pymol_Visualize_Path(graph=intersection_graph,
                                                                                             pdb_file=self.pdb)

                            # ----------------------> 3D NETWORK VISUALIZATION USING PYMOL / START <---------------------- #
                            Protein3DNetworkView.show_energy_dissipation(response_time_file_path=self.retime_file)
                            for arrow_coord in arrows_cordinates:
                                Protein3DNetworkView.create_directed_arrows(atom1=arrow_coord[0], atom2=arrow_coord[1],
                                                                            radius=0.05,
                                                                            gap=0.4, hradius=0.4, hlength=0.8,
                                                                            color='green')
                            for node in intersection_node_list:
                                resID_of_node = int(''.join(list(filter(str.isdigit, node))))
                                Protein3DNetworkView.resi_label_add('resi ' + str(resID_of_node))

                            # MAKE PYMOL VISUALIZATION BETTER
                            Protein3DNetworkView._pymol.cmd.set('cartoon_oval_length', 0.8)  # default is 1.20)
                            Protein3DNetworkView._pymol.cmd.set('cartoon_oval_width', 0.2)
                            Protein3DNetworkView._pymol.cmd.center(selection="all", state=0, origin=1, animate=0)
                            Protein3DNetworkView._pymol.cmd.zoom('all', buffer=0.0, state=0, complete=0)
                            Protein3DNetworkView.update()
                            Protein3DNetworkView.show()

                            # ----------------------> 2D NETWORK VISUALIZATION USING visJS / START <---------------------- #
                            # self.load_nx_to_VisJS_2D_Network(intersection_gml_file=intersection_graph)

                        except Exception as error:
                            print("Problem: ", error)
                            exc_type, exc_obj, exc_tb = sys.exc_info()
                            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
                            print(exc_type, fname, exc_tb.tb_lineno)

                        all = ''
                        for j in self.log_holder:
                            if j != '':
                                all = all + '\n' + j
                                print("all: ", all)
                        Message_Boxes.Information_message(self, "DONE !", all, Style.MessageBox_stylesheet)
                        del self.log_holder, self.network_holder

                    # elif len(intersection_graph.nodes()) == 0:
                    #     print("NO INTERSECTION GRAPH")

                    else:
                        Message_Boxes.Information_message(self, "DONE !", "There is no Intersection Network :(",
                                                          Style.MessageBox_stylesheet)
            except Exception as e:
                print("Problem: ", e)
                exc_type, exc_obj, exc_tb = sys.exc_info()
                fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
                print(exc_type, fname, exc_tb.tb_lineno)
        else:
            print("There is no suitable Graph for your search parameters")

        #
        # finally:
        #     # CLOSE THE ALREADY OPENED POOL
        #     pool.close()
        #     pool.join()

    def show_shortest_paths_on_3D_ProteinView(self, item, PyMOL_Widget):
        processed_path = [x.strip() for x in item.text().split('-->')]
        shortest_path_arrow_coords = Shortest_Path_Visualize(pdb_file=self.pdb, selected_path=processed_path)

        for arrow_coord in shortest_path_arrow_coords:
            PyMOL_Widget.create_directed_arrows(atom1=arrow_coord[0], atom2=arrow_coord[1],
                                                radius=0.05,
                                                gap=0.4, hradius=0.4, hlength=0.8,
                                                color='green', shortest_path=True)
        for node in processed_path:
            resID_of_node = int(''.join(list(filter(str.isdigit, node))))
            PyMOL_Widget.resi_label_add('resi ' + str(resID_of_node))

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
            options = QFileDialog.Options()
            options |= QFileDialog.DontUseNativeDialog
            self.response_filename, _ = QFileDialog.getOpenFileName(self, "Select The *.csv File", str(os.getcwd()),
                                                                    "Response_Time_Files (*.csv)", str(options))

            self.response_time_file_path = self.response_filename
            self.response_filename = os.path.splitext(os.path.basename(self.response_filename))

            if self.response_filename[1] == '.csv':
                self.response_time_lineEdit.setText(self.response_time_file_path)
                return True, self.response_time_file_path

            elif self.response_filename[1] != "":
                Message_Boxes.Critical_message(self, "Error", "this is not a valid response time file",
                                               Style.MessageBox_stylesheet)

        except Exception as exp:
            Message_Boxes.Warning_message(self, "Fatal Error!", str(exp), Style.MessageBox_stylesheet)

    def browse_bound_form_pdbFile(self):
        """
            The function provides Main GUI / Upload button activity for select pdb file indicated by the user
        """
        try:
            options = QFileDialog.Options()
            options |= QFileDialog.DontUseNativeDialog
            self.boundForm_pdb_filename, _ = QFileDialog.getOpenFileName(self, "Show The *pdb File", str(os.getcwd()),
                                                                         "pdb Files (*.pdb)", str(options))

            self.boundForm_pdb_path = self.boundForm_pdb_filename
            self.boundForm_pdb_filename = os.path.splitext(os.path.basename(self.boundForm_pdb_filename))

            if self.boundForm_pdb_filename[1] == '.pdb':
                self.boundForm_pdb_lineedit.setText(self.boundForm_pdb_path)
                return True, self.boundForm_pdb_path

            elif self.boundForm_pdb_filename[1] != "":
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

    # ########################################### ANALYSIS WINDOW FUNCTIONS ############################################

    # ######################################### PERTURBATION WINDOW FUNCTIONS ##########################################
    def maximum_thread_of_system(self):
        self.Number_CPU_spinBox.setMaximum(mp.cpu_count())
        self.Number_of_thread_for_network_spinBox.setMaximum(mp.cpu_count())

    def precision_combobox_settings(self, eq_md_indexes, per_md_indexes):

        if eq_md_indexes is not None:
            if str(self.equ_platform_comboBox.currentText()) == 'CPU':
                self.eq_precision_comboBox.model().item(eq_md_indexes['single']).setEnabled(False)
                self.eq_precision_comboBox.model().item(eq_md_indexes['mixed']).setEnabled(True)
                self.eq_precision_comboBox.model().item(eq_md_indexes['double']).setEnabled(False)
                self.eq_precision_comboBox.setCurrentIndex(eq_md_indexes['mixed'])

            if str(self.equ_platform_comboBox.currentText()) == 'Reference':
                self.eq_precision_comboBox.model().item(eq_md_indexes['single']).setEnabled(False)
                self.eq_precision_comboBox.model().item(eq_md_indexes['mixed']).setEnabled(False)
                self.eq_precision_comboBox.model().item(eq_md_indexes['double']).setEnabled(True)
                self.eq_precision_comboBox.setCurrentIndex(eq_md_indexes['double'])

            if str(self.equ_platform_comboBox.currentText()) in ['CUDA', 'OpenCL']:
                self.eq_precision_comboBox.model().item(eq_md_indexes['single']).setEnabled(True)
                self.eq_precision_comboBox.model().item(eq_md_indexes['mixed']).setEnabled(True)
                self.eq_precision_comboBox.model().item(eq_md_indexes['double']).setEnabled(True)
                self.eq_precision_comboBox.setCurrentIndex(eq_md_indexes['single'])

        if per_md_indexes is not None:
            if str(self.per_platform_comboBox.currentText()) == 'CPU':
                self.per_precision_comboBox.model().item(per_md_indexes['single']).setEnabled(False)
                self.per_precision_comboBox.model().item(per_md_indexes['mixed']).setEnabled(True)
                self.per_precision_comboBox.model().item(per_md_indexes['double']).setEnabled(False)
                self.per_precision_comboBox.setCurrentIndex(per_md_indexes['mixed'])

            if str(self.per_platform_comboBox.currentText()) == 'Reference':
                self.per_precision_comboBox.model().item(per_md_indexes['single']).setEnabled(False)
                self.per_precision_comboBox.model().item(per_md_indexes['mixed']).setEnabled(False)
                self.per_precision_comboBox.model().item(per_md_indexes['double']).setEnabled(True)
                self.per_precision_comboBox.setCurrentIndex(per_md_indexes['double'])

            if str(self.per_platform_comboBox.currentText()) in ['CUDA', 'OpenCL']:
                self.per_precision_comboBox.model().item(per_md_indexes['single']).setEnabled(True)
                self.per_precision_comboBox.model().item(per_md_indexes['mixed']).setEnabled(True)
                self.per_precision_comboBox.model().item(per_md_indexes['double']).setEnabled(True)
                self.per_precision_comboBox.setCurrentIndex(per_md_indexes['single'])

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
        from pathlib import Path
        current_path = os.path.dirname(os.path.realpath(__file__))
        path = os.path.join(Path(current_path).parent, 'Download', '4htm.pdb')
        output_directory = os.path.join(Path(current_path).parent, 'output')

        if os.path.exists(path):
            self.upload_pdb_lineEdit.setText(path)
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
                self.selected_residues_listWidget.addItem('SER345')
            self.R_factor_lineEdit.setText("3, 4")

            self.run_duration_spinBox.blockSignals(True)
            self.run_duration_doubleSpinBox.blockSignals(True)
            self.Number_of_steps_spinBox.blockSignals(True)
            self.run_duration_spinBox.setValue(int(1000))
            self.run_duration_doubleSpinBox.setValue(float(1.000))
            self.Number_of_steps_spinBox.setValue(int(500000))
            self.run_duration_spinBox.blockSignals(False)
            self.run_duration_doubleSpinBox.blockSignals(False)
            self.Number_of_steps_spinBox.blockSignals(False)

    def export_workspace(self):
        try:
            # default_dir = os.path.join(os.path.join(os.path.expanduser('~')), 'Desktop')

            options = QFileDialog.Options()
            options |= QFileDialog.DontUseNativeDialog
            output_file_directory, _ = QFileDialog.getSaveFileName(self, "Save configuration file", "",
                                                                   "Config File Name (*.yaml)", options=options)
            template = config_template()

            # INPUTS
            template['Inputs']['topology'] = self.upload_pdb_lineEdit.text().strip()
            template['Inputs']['output directory'] = self.Output_Folder_textEdit.toPlainText().strip()
            template['Inputs']['run duration'] = float(self.run_duration_doubleSpinBox.value())
            template['Inputs']['run duration unit'] = self.long_simulation_time_unit.currentText()
            template['Inputs']['threading number'] = int(self.Number_CPU_spinBox_2.value())
            template['Inputs']['all cpu is active'] = self.All_CPU_checkBox.isChecked()
            template['Inputs']['r factor'] = self.R_factor_lineEdit.text().strip()
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
            template['Simulation']['equilibrium time step'] = self.integrator_time_step.toPlainText().strip()
            template['Simulation']['equilibrium time step unit'] = self.integrator_time_step_unit.currentText()
            template['Simulation'][
                'equilibrium additional integrator parameters'] = self.Additional_Integrators_checkBox.isChecked()
            template['Simulation']['friction'] = self.friction_textEdit.toPlainText().strip()
            template['Simulation']['temperature'] = self.temperature_textEdit.toPlainText().strip()
            template['Simulation']['nonbonded method'] = self.nonBounded_Method_comboBox.currentText()
            template['Simulation']['constraints'] = self.system_constraints_comboBox.currentText()
            template['Simulation']['rigid water is active'] = self.rigid_water_checkBox.isChecked()
            template['Simulation']['nonbonded cutoff'] = self.nonbounded_CutOff_textEdit.toPlainText()
            template['Simulation']['switching distance'] = self.switching_distance_textEdit.toPlainText()
            template['Simulation']['use switching distance'] = self.use_switching_checkBox.isChecked()
            template['Simulation']['number of simulation steps'] = self.Number_of_steps_spinBox.value()
            template['Simulation']['minimize'] = self.minimize_checkBox.isChecked()
            template['Simulation']['minimization max iterations'] = self.Max_minimize_iter_textEdit.toPlainText()
            template['Simulation']['Equilibrate'] = self.equilubrate_checkBox.isChecked()
            template['Simulation']['Equilibration steps'] = self.Max_equilubrate_steps_textEdit.toPlainText()
            template['Simulation']['DCD report'] = self.DCD_Reporter_checkBox.isChecked()
            template['Simulation']['XTC report'] = self.XTC_Reporter_checkBox.isChecked()
            template['Simulation']['state data report'] = self.State_Data_Reporter_checkBox.isChecked()
            template['Simulation']['DCD writing frequency'] = self.DCD_write_freq_textEdit.toPlainText()
            template['Simulation']['DCD output name'] = self.DCD_Output_Name_textEdit.toPlainText()
            template['Simulation']['XTC writing frequency'] = self.XTC_write_freq_textEdit.toPlainText()
            template['Simulation']['XTC output name'] = self.XTC_Output_Name_textEdit.toPlainText()
            template['Simulation']['state data frequency'] = self.StateData_frequency_textEdit.toPlainText()

            write_output_configuration_file(file_path=output_file_directory, template_yml=template)

        except Exception as err:
            print("ERR: ", err)
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
            self.R_factor_lineEdit.setText(str(template['Inputs']['r factor']))
            self.selected_residues_listWidget.addItems(template['Inputs']['selected residue'])
            self.run_duration_spinBox.setValue(int(template['Inputs']['perturbation run duration']))
            self.perturb_time_unit_comboBox.setCurrentIndex(
                self.perturb_time_unit_comboBox.findText(
                    str(template['Inputs']['perturbation run duration unit']), QtCore.Qt.MatchFixedString))
            self.Number_CPU_spinBox.setValue(int(template['Inputs']['perturbation threading number']))
            self.perturbation_All_CPU_checkBox.setChecked(template['Inputs']['perturbation all cpu is active'])

        except Exception as Err:
            print(Err)
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
                                                                    self.Output_Folder_textEdit.toPlainText(), ph=7,
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
                                                                     ph=7, chains_to_remove=None)

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

    def Stochastic_changed(self):
        """
             The function provides enable or disable of Stochastic Parameter menu according to İntegrator kind
        """
        self.kind_of_integrator = self.integrator_kind_comboBox.currentText()

        if self.kind_of_integrator in ["Langevin", "Brownian"]:
            self.Additional_Integrator_groupBox.setEnabled(True)
            self.Additional_Integrators_checkBox.setChecked(True)
            self.Additional_Integrators_checkBox.setEnabled(True)
        else:
            self.Additional_Integrator_groupBox.setEnabled(False)
            self.Additional_Integrators_checkBox.setChecked(False)
            self.Additional_Integrators_checkBox.setEnabled(False)

    @staticmethod
    def minimize_Step_isVisible(self):
        if not self.minimize_checkBox.isChecked():
            self.Max_minimize_iter_textEdit.setEnabled(False)
        else:
            self.Max_minimize_iter_textEdit.setEnabled(True)

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
        current_integrator_time_step_value = float(self.integrator_time_step.toPlainText())

        if current_time_unit == 'nanosecond':
            new_step = int((current_step / current_integrator_time_step_value) * 1000000)

        if current_time_unit == 'picosecond':
            new_step = int((current_step / current_integrator_time_step_value) * 1000)

        self.Number_of_steps_spinBox.setValue(new_step)

    def number_of_steps_changed_from_advanced(self):
        global new_time
        current_step = int(self.Number_of_steps_spinBox.value())  # 1 ns
        current_time_unit = self.long_simulation_time_unit.currentText()  # ns
        current_integrator_time_step_value = float(self.integrator_time_step.toPlainText())  # 2 fs

        if current_time_unit == 'nanosecond':
            new_time = float((current_step * current_integrator_time_step_value) / 1000000)

        if current_time_unit == 'picosecond':
            new_time = float((current_step * current_integrator_time_step_value) / 1000)

        self.run_duration_doubleSpinBox.setValue(new_time)


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
            :param file_pathway: pathway for manipulating your fetched pdb files
            :param chains_to_remove: Selected chains will be deleted
            :param ph: Selected pH value will be apply to the structure's Hydrogens
        Returns:
            :param output_path: the manipulated pdb file will return as full path if specified
                                otherwise will return already exist path
        """
        ## get name of pdb file ##
        name_of_pdb = os.path.basename(file_pathway).split('.')[0]

        print("Creating PDBFixer...")
        fixer = PDBFixer(file_pathway)
        print("Finding missing residues...")

        if chains_to_remove is not None:
            print("toDelete: %s" % chains_to_remove)
            fixer.removeChains(chainIds=chains_to_remove)

        fixer.findMissingResidues()

        chains = list(fixer.topology.chains())
        keys = fixer.missingResidues.keys()
        for key in list(keys):
            chain = chains[key[0]]
            if key[1] == 0 or key[1] == len(list(chain.residues())):
                print("ok")
                del fixer.missingResidues[key]

        print("Finding nonstandard residues...")
        fixer.findNonstandardResidues()
        print("Replacing nonstandard residues...")
        fixer.replaceNonstandardResidues()
        print("Removing heterogens...")
        fixer.removeHeterogens(keepWater=False)
        """
        print("Finding missing atoms...")
        fixer.findMissingAtoms()
        print("Adding missing atoms...")
        fixer.addMissingAtoms()
        print("Adding missing hydrogens...")
        fixer.addMissingHydrogens(pH=ph)
        """
        print("Writing PDB file...")

        ####  FOR DELETE WITH MODELLER USE FOLLOWING SCRIPT

        # if chains_to_remove is not None:
        #     toDelete = [r for r in modeller.topology.chains() if r.id in chains_to_remove]
        #     modeller.delete(toDelete)

        if output_path != "":
            PDBFile.writeFile(
                fixer.topology,
                fixer.positions,
                open(os.path.join(output_path, "%s_fixed_ph%s.pdb" % (name_of_pdb, ph)),
                     "w"),
                keepIds=True)
            return os.path.join(output_path, "%s_fixed_ph%s.pdb" % (name_of_pdb, ph))

        if output_path == "":
            new_outpath_dir = os.path.dirname(file_pathway)
            new_outpath = os.path.join(new_outpath_dir, "%s_fixed_ph%s.pdb" % (name_of_pdb, ph))
            PDBFile.writeFile(
                fixer.topology,
                fixer.positions,
                open(new_outpath, "w"), keepIds=True)

            return new_outpath

        # # Remove the ligand and write a pdb file
        # fixer.removeHeterogens(True)
        # PDBFile.writeFile(
        #     fixer.topology,
        #     fixer.positions,
        #     open(os.path.join('/home/enaz/Desktop/MDPERTOOL_v01/Output',
        #                       "%s_fixed_ph%s_apo.pdb" % ('pdbid', 7)), "w"),
        #     keepIds=True)

#
# pdb = 'C:\\Users\\HIbrahim\\Desktop\\MDPERTOOL_v01\\Download\\1gg1.pdb'
# out = 'C:\\Users\\HIbrahim\\Desktop\\MDPERTOOL_v01\\Download'
# pdb_Tools().fetched_pdb_fix(pdb, output_path=out, ph=7, chains_to_remove=['B', 'C'])
