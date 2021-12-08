import os
from PySide2.QtWidgets import QFileDialog, QMessageBox



class browse_file:

    def __init__(self):
        self.upload_pdb_textEdit = None

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
            self.upload_pdb_textEdit.setText(self.pdb_path)
            return True, self.pdb_path
        elif self.pdb_filename[1] != "":
            QMessageBox.critical(self, "Error", "this is not a pdb file")

    def PDB_ID_textEdit(self):
        """
             The function provides enable or disable of Upload Button according to the entry of user
        """
        self.text_edit = self.PDB_ID_textEdit.toPlainText()
        if self.text_edit == "":
            self.upload_pdb_Button.setEnabled(True)
            self.upload_pdb_textEdit.setEnabled(True)
            self.label_32.setEnabled(True)
        else:
            self.upload_pdb_Button.setEnabled(False)
            self.upload_pdb_textEdit.setEnabled(False)
            self.label_32.setEnabled(False)

    def Stochastic_changed(self):
        """
             The function provides enable or disable of Stochastic Parameter menu according to İntegrator kind
        """
        self.kind_of_integrator = self.integrator_kind_comboBox.currentText()
        if (self.kind_of_integrator in ["Langevin", "Brownian"]):
            self.Additional_Integrator_groupBox.setEnabled(False)
            self.Additional_Integrators_checkBox.setChecked(False)
            self.Additional_Integrators_checkBox.setEnabled(False)
        else:
            self.Additional_Integrator_groupBox.setEnabled(True)
            self.Additional_Integrators_checkBox.setChecked(True)
            self.Additional_Integrators_checkBox.setEnabled(True)

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

    def State_Data_Reporter_Changed(self):
        if not self.State_Data_Reporter_checkBox.isChecked():
            self.State_Data_Reporter_Options_groupBox.setEnabled(False)

        if self.State_Data_Reporter_checkBox.isChecked():
            self.State_Data_Reporter_Options_groupBox.setEnabled(True)

    def output_file(self):
        try:
            options = QFileDialog.Options()
            options |= QFileDialog.DontUseNativeDialog
            self.output_file = QFileDialog.getExistingDirectory(options=options)
            self.Output_Folder_textEdit.setText(self.output_file)
            return True
        except:
            return False

    def available_platforms(self):
        """
            :return: The function will return available platforms on your system for OpenMM Engine
        """
        import simtk.openmm
        return [simtk.openmm.Platform.getPlatform(index).getName() for index in
                range(simtk.openmm.Platform.getNumPlatforms())]

    def fill_residue_combobox(self, pdb_path):
        from prody.proteins.pdbfile import parsePDB
        self.combobox = []
        self.pdb_path = pdb_path
        self.pdb = parsePDB(self.pdb_path)
        self.protein = self.pdb.select('protein')
        for model in self.protein.getHierView():
            for chain in model:
                # print(chain)
                # self.combobox.append(str(model).split(" ")[1] + str(chain))
                self.combobox.append(str(chain).replace(" ", "") + str(model).split(" ")[1])

                # self.combobox_2.append(str(model).split(" ")[1] + str(chain))
        for i in self.combobox:
            self.res1_comboBox.addItem(str(i))
            self.res2_comboBox.addItem(str(i))
        self.res1_comboBox.clear()  # delete all items from comboBox
        self.res1_comboBox.addItems(self.combobox)  # add the actual content of self.comboData
        self.res2_comboBox.clear()  # delete all items from comboBox
        self.res2_comboBox.addItems(self.combobox)  # add the actual content of self.comboData
