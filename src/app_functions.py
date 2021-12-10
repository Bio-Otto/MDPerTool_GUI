from ui_main import *
import networkx as nx
from PySide2.QtWidgets import QFileDialog, QDialog, QTableWidgetItem
from PySide2.QtCore import Qt, Slot
import csv
from pathlib import Path
from pdbfixer import PDBFixer
from os import path
from urllib.request import urlretrieve
from openmm.app import *
from src.checkBox_menu import *
from src.message import Message_Boxes
import multiprocessing as mp
from analysis.pdbsum_conservation_puller import get_conservation_scores
from analysis.createRNetwork import (Multi_Task_Engine, intersection_of_directed_networks, Pymol_Visualize_Path)


class Helper_Functions():

    def fill_residue_combobox(self, pdb_path):
        print(pdb_path)
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
            self.plot_signal.plot_network.emit()
            self.progress.cancel()

    def print_output(self, s):
        self.network_holder.append(s[0])
        self.log_holder.append(s[1])

    def calculate_intersection_network(self):
        global intersection_graph, output_folder_directory, network_holder

        # try:
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

        if self.node_threshold_use_condition:
            self.node_threshold = None

        verbose_condition = True  # True or False
        target_residues = [self.selected_target_residues_listWidget.item(x).text()[:-1]
                           for x in range(self.selected_target_residues_listWidget.count())]  # None or residue list

        use_conservation = self.use_conservation_checkBox.isChecked()

        pdb_id = self.conservation_PDB_ID_lineEdit.text()  # FREE OR BOUND FORM OF PDB, like '2EB8'
        chain = self.conservation_pdb_chain_id_lineedit.text()  # FOR PULLING CONS. SCORES INDICATE CHAIN ID OF PDB
        conservation_threshold = self.conserv_score_doubleSpinBox.value()
        save_conservation_scores = False

        visualize_on_PyMol = True  # Networkx Graph Visulization on PyMol
        visualize_on_VisJS = True  # Networkx Graph Visulization on VisJS
        self.create_output = True  # Supports True or False Conditions for creation of all networks (*.gml) on a folder
        just_visualize = False  # If you have already calculated network you can directly visualize it.

        general_output_folder = os.path.join(self.output_directory, 'network_outputs')
        Path(general_output_folder).mkdir(parents=True, exist_ok=True)

        folder_name = "output_%s" % self.source
        output_folder_directory = os.path.join(general_output_folder, folder_name)
        Path(output_folder_directory).mkdir(parents=True, exist_ok=True)

        # pool = mp.Pool(number_of_threads)
        engine = Multi_Task_Engine(pdb_file=self.pdb, cutoff=self.cutoff, reTimeFile=self.retime_file,
                                   source=self.source,
                                   node_threshold=self.node_threshold, verbose=verbose_condition,
                                   outputFileName=self.outputFileName, write_outputs=self.create_output,
                                   output_directory=output_folder_directory)

        engine.calculate_general_network()
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

        # time.sleep(0.1)
        # net, log = zip(*pool.map(engine, target_residues))  # ########################################

    def plot_networks(self):
        print(len(self.network_holder))

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
            intersection_graph = intersection_of_directed_networks(clean_graph_list)
            if self.create_output:
                nx.write_gml(intersection_graph, os.path.join(output_folder_directory, 'intersection_graph.gml'))

        else:
            print("There is no suitable Graph for your search parameters")

        try:
            if len(intersection_graph.nodes()) > 0:
                try:

                    arrows_cordinates, intersection_node_list = Pymol_Visualize_Path(graph=intersection_graph,
                                                                                     pdb_file=self.pdb)

                    # ----------------------> 3D NETWORK VISUALIZATION USING PYMOL / START <---------------------- #
                    self.Protein3DNetworkView.show_energy_dissipation(response_time_file_path=self.retime_file)
                    for arrow_coord in arrows_cordinates:
                        self.Protein3DNetworkView.create_directed_arrows(atom1=arrow_coord[0], atom2=arrow_coord[1],
                                                                         radius=0.05,
                                                                         gap=0.4, hradius=0.4, hlength=0.8, color='green')
                    for node in intersection_node_list:
                        resID_of_node = int(''.join(list(filter(str.isdigit, node))))
                        self.Protein3DNetworkView.resi_label_add('resi ' + str(resID_of_node))

                    # MAKE PYMOL VISUALIZATION BETTER
                    self.Protein3DNetworkView._pymol.cmd.set('cartoon_oval_length', 0.8)  # default is 1.20)
                    self.Protein3DNetworkView._pymol.cmd.set('cartoon_oval_width', 0.2)
                    self.Protein3DNetworkView._pymol.cmd.center(selection="all", state=0, origin=1, animate=0)
                    self.Protein3DNetworkView._pymol.cmd.zoom('all', buffer=0.0, state=0, complete=0)
                    self.Protein3DNetworkView.update()
                    self.Protein3DNetworkView.show()

                    # ----------------------> 2D NETWORK VISUALIZATION USING visJS / START <---------------------- #
                    self.load_nx_to_VisJS_2D_Network(intersection_gml_file=intersection_graph)

                except Exception as error:
                    print("Problem: ", error)
                    exc_type, exc_obj, exc_tb = sys.exc_info()
                    fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
                    print(exc_type, fname, exc_tb.tb_lineno)

                all = ''
                for j in self.log_holder:
                    all = all + '\n' + j
                Message_Boxes.Information_message(self, "DONE !", all, Style.MessageBox_stylesheet)
                del self.log_holder, self.network_holder
            else:
                Message_Boxes.Information_message(self, "DONE !", "There is no Intersection Network :(", Style.MessageBox_stylesheet)
        except Exception as e:
            print("Problem: ", e)
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
            print(exc_type, fname, exc_tb.tb_lineno)
        #
        # finally:
        #     # CLOSE THE ALREADY OPENED POOL
        #     pool.close()
        #     pool.join()

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

    # ########################################### ANALYSIS WINDOW FUNCTIONS ############################################

    # ######################################### PERTURBATION WINDOW FUNCTIONS ##########################################
    def maximum_thread_of_system(self):
        self.Number_CPU_spinBox.setMaximum(mp.cpu_count())
        self.Number_of_thread_for_network_spinBox.setMaximum(mp.cpu_count())

    def output_file(self):
        try:
            options = QFileDialog.Options()
            options |= QFileDialog.DontUseNativeDialog
            output_file = QFileDialog.getExistingDirectory(options=options)
            self.Output_Folder_textEdit.setText(output_file)
            return True

        except Exception as ins:
            print(ins)
            return False

    def browse_pdbFile(self):
        """
            The function provides Main GUI / Upload button activity for select pdb file indicated by the user
        """
        try:
            options = QFileDialog.Options()
            options |= QFileDialog.DontUseNativeDialog
            self.pdb_filename, _ = QFileDialog.getOpenFileName(self, "Show The *pdb File", str(os.getcwd()),
                                                               "pdb Files (*.pdb)", str(options))

            self.pdb_path = self.pdb_filename
            self.pdb_filename = os.path.splitext(os.path.basename(self.pdb_filename))

            if self.pdb_filename[1] == '.pdb':
                return True, self.pdb_path

            elif self.pdb_filename[1] != "":
                Message_Boxes.Critical_message(self, "Error", "this is not a pdb file", Style.MessageBox_stylesheet)

        except Exception as exp:
            Message_Boxes.Warning_message(self, "Fatal Error!", str(exp), Style.MessageBox_stylesheet)

    @staticmethod
    def PDB_ID_lineEdit(self):
        """
             The function provides enable or disable of Upload Button according to the entry of user
        """
        self.text_edit = self.PDB_ID_lineEdit.text()
        if self.text_edit == "":
            self.upload_pdb_Button.setEnabled(True)
            self.upload_pdb_textEdit.setEnabled(True)
            self.label_32.setEnabled(True)
        else:
            self.upload_pdb_Button.setEnabled(False)
            self.upload_pdb_textEdit.setEnabled(False)
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
                            print(delete_chains)
                            print("fetched pdb: %s" % fetched_pdb)

                            self.upload_pdb_textEdit.setText(fetched_pdb)
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

                            self.upload_pdb_textEdit.setText(modified_pdb)

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
        self.platform_list_on_the_program = [self.platform_comboBox.itemText(i) for i in
                                             range(self.platform_comboBox.count())]

        for item_no, i in enumerate(self.platform_list_on_the_program):
            if i not in self.platforms:
                print(item_no)
                self.platform_comboBox.model().item(int(item_no)).setEnabled(False)
                self.platform_comboBox.setCurrentIndex(item_no + 1)

            if i in self.platforms:
                self.platform_comboBox.setItemData(item_no, str("Estimated Speed For This Devices Is "
                                                                + str(self.plt_speeds[i])), Qt.ToolTipRole)

    @staticmethod
    def platform_comboBox_Changed(self):
        if self.platform_comboBox.currentText() in ["CPU", "Reference"]:
            self.Device_Number_comboBox.setEnabled(False)
            self.Device_ID_checkBox.setEnabled(False)

        else:
            self.Device_Number_comboBox.setEnabled(True)
            self.Device_ID_checkBox.setEnabled(True)

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
            print(Download_folder)

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
