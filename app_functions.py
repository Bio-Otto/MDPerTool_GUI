from PyQt5.QtWidgets import QFileDialog, QWidget
from PyQt5.QtWidgets import QMessageBox

# from PyQt5.QtWidgets import *
# from PyQt5 import QtCore, QtGui, QtWidgets
# from PyQt5.QtCore import (QCoreApplication, QPropertyAnimation, QDate, QDateTime, QMetaObject, QObject, QPoint, QRect,
#                           QSize, QTime, QUrl, Qt, QEvent)
# from PyQt5.QtGui import (QBrush, QColor, QConicalGradient, QCursor, QFont, QFontDatabase, QIcon, QKeySequence,
#                          QLinearGradient, QPalette, QPainter, QPixmap, QRadialGradient)
#

import gzip

from PyQt5 import QtCore
from PyQt5.QtCore import pyqtSlot
from PyQt5.QtWidgets import QMessageBox
from pdbfixer import PDBFixer
from simtk.openmm import *
from simtk.openmm.app import *

from checkBox_menu import *
from os import path
from urllib.request import urlretrieve
from checkBox_menu import *
from ui_main import *


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
        import simtk.openmm
        return [simtk.openmm.Platform.getPlatform(index).getName() for index in
                range(simtk.openmm.Platform.getNumPlatforms())]


class Functions(MainWindow):

    def output_file(self):
        try:
            options = QFileDialog.Options()
            options |= QFileDialog.DontUseNativeDialog
            output_file = QFileDialog.getExistingDirectory(options=options)
            self.Output_Folder_textEdit.setText(output_file)
            return True
        except:
            return False

    def browse_pdbFile(self):
        """
            The function provides Main GUI / Upload button activity for select pdb file indicated by the user
        """

        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        self.pdb_filename, _ = QFileDialog.getOpenFileName(self, "Show The *pdb File", "", "pdb Files (*.pdb)",
                                                           str(options))
        self.pdb_path = self.pdb_filename

        self.pdb_filename = os.path.basename(self.pdb_filename)
        self.pdb_filename = os.path.splitext(self.pdb_filename)

        if self.pdb_filename[1] == '.pdb':
            # self.upload_pdb_textEdit.setText(self.pdb_path)
            return True, self.pdb_path
        elif self.pdb_filename[1] != "":
            QMessageBox.critical(self, "Error", "this is not a pdb file")

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

    def Fetch_PDB_File(self):
        global selected_chains
        try:
            fetch_pdb_ID = self.PDB_ID_lineEdit.text()
            if len(fetch_pdb_ID) == 4:

                fetch_result = pdb_Tools.fetch_pdb(self, fetch_pdb_ID)

                if path.exists(fetch_result):
                    fixer = PDBFixer(fetch_result)
                    fixer.removeHeterogens(keepWater=False)

                    modeller = Modeller(fixer.topology, fixer.positions)
                    chains = [r.id for r in modeller.topology.chains()]

                    checked_list = ChecklistDialog('Select the chain (s) to be used in the system', chains,
                                                   checked=True)
                    if checked_list.exec_() == QtWidgets.QDialog.Accepted:
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
                        self.res2_comboBox.addItem(str(i))
                    self.res1_comboBox.clear()  # delete all items from comboBox
                    self.res1_comboBox.addItems(self.combobox)  # add the actual content of self.comboData
                    self.res2_comboBox.clear()  # delete all items from comboBox
                    self.res2_comboBox.addItems(self.combobox)  # add the actual content of self.comboData
                    InputFile(fetch_result)

            if len(fetch_pdb_ID) != 4:
                QMessageBox.critical(self, 'Wrong pdb id', 'PDB ID should be provided as 4 letters')

        except Exception as instance:
            QMessageBox.critical(self, 'An error occurred while fetching the pdb file.', repr(instance))

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

    def minimize_Step_isVisible(self):
        if not self.minimize_checkBox.isChecked():
            self.Max_minimize_iter_textEdit.setEnabled(False)
        else:
            self.Max_minimize_iter_textEdit.setEnabled(True)

    def DCD_Reporter_Changed(self):
        if not self.DCD_Reporter_checkBox.isChecked():
            self.DCD_Reporter_Options_groupBox.setEnabled(False)

        if self.DCD_Reporter_checkBox.isChecked():
            self.DCD_Reporter_Options_groupBox.setEnabled(True)

    def XTC_Reporter_Changed(self):
        if not self.XTC_Reporter_checkBox.isChecked():
            self.XTC_Reporter_Options_groupBox.setEnabled(False)

        if self.XTC_Reporter_checkBox.isChecked():
            self.XTC_Reporter_Options_groupBox.setEnabled(True)

    def State_Data_Reporter_Changed(self):
        if not self.State_Data_Reporter_checkBox.isChecked():
            self.State_Data_Reporter_Options_groupBox.setEnabled(False)

        if self.State_Data_Reporter_checkBox.isChecked():
            self.State_Data_Reporter_Options_groupBox.setEnabled(True)

    def Send_Available_Platforms_to_GUI(self):
        self.platforms = Helper_Functions.available_platforms(self)
        self.platform_list_on_the_program = [self.platform_comboBox.itemText(i) for i in
                                             range(self.platform_comboBox.count())]

        for item_no, i in enumerate(self.platform_list_on_the_program):
            if i not in self.platforms:
                print(item_no)
                self.platform_comboBox.model().item(int(item_no)).setEnabled(False)
                self.platform_comboBox.setCurrentIndex(item_no + 1)

    def platform_comboBox_Changed(self):
        if self.platform_comboBox.currentText() in ["CPU", "Reference"]:
            self.Device_Number_comboBox.setEnabled(False)
            self.Device_ID_checkBox.setEnabled(False)

        else:
            self.Device_Number_comboBox.setEnabled(True)
            self.Device_ID_checkBox.setEnabled(True)


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
            QMessageBox.critical(self, 'An error occurred while fetching the pdb file.', repr(instance))
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
            return output_path

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
