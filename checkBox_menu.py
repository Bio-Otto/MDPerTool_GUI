from PyQt5 import QtWidgets, QtGui, QtCore
import sys
from pdbfixer import PDBFixer
from simtk.openmm import *
from simtk.openmm.app import *


class ChecklistDialog(QtWidgets.QDialog):

    def __init__(self,
                 name,
                 stringlist=None,
                 checked=False,
                 icon=None,
                 parent=None,
                 ):
        super(ChecklistDialog, self).__init__(parent)

        self.name = name
        self.icon = icon
        self.model = QtGui.QStandardItemModel()
        self.listView = QtWidgets.QListView()

        for string in stringlist:
            item = QtGui.QStandardItem(string)
            item.setCheckable(True)
            check = \
                (QtCore.Qt.Checked if checked else QtCore.Qt.Unchecked)
            item.setCheckState(check)
            self.model.appendRow(item)

        self.listView.setModel(self.model)

        self.okButton = QtWidgets.QPushButton('OK')
        self.cancelButton = QtWidgets.QPushButton('Cancel')
        self.selectButton = QtWidgets.QPushButton('Select All')
        self.unselectButton = QtWidgets.QPushButton('Unselect All')

        hbox = QtWidgets.QHBoxLayout()
        hbox.addStretch(1)
        hbox.addWidget(self.okButton)
        hbox.addWidget(self.cancelButton)
        hbox.addWidget(self.selectButton)
        hbox.addWidget(self.unselectButton)

        vbox = QtWidgets.QVBoxLayout(self)
        vbox.addWidget(self.listView)
        vbox.addStretch(1)
        vbox.addLayout(hbox)

        self.setWindowTitle(self.name)
        if self.icon:
            self.setWindowIcon(self.icon)

        self.okButton.clicked.connect(self.onAccepted)
        self.cancelButton.clicked.connect(self.reject)
        self.selectButton.clicked.connect(self.select)
        self.unselectButton.clicked.connect(self.unselect)

    def onAccepted(self):
        self.choices = [self.model.item(i).text() for i in
                        range(self.model.rowCount())
                        if self.model.item(i).checkState()
                        == QtCore.Qt.Checked]
        self.accept()

    def select(self):
        for i in range(self.model.rowCount()):
            item = self.model.item(i)
            item.setCheckState(QtCore.Qt.Checked)

    def unselect(self):
        for i in range(self.model.rowCount()):
            item = self.model.item(i)
            item.setCheckState(QtCore.Qt.Unchecked)




# fruits = [
#     'Banana',
#     'Apple',
#     'Elderberry',
#     'Clementine',
#     'Fig',
#     'Guava',
#     'Mango',
#     'Honeydew Melon',
#     'Date',
#     'Watermelon',
#     'Tangerine',
#     'Ugli Fruit',
#     'Juniperberry',
#     'Kiwi',
#     'Lemon',
#     'Nectarine',
#     'Plum',
#     'Raspberry',
#     'Strawberry',
#     'Orange',
# ]

# fetch_result = '/home/enaz/Desktop/MDPERTOOL_v01/Download/7adb.pdb'
# #
# fixer = PDBFixer(fetch_result)
# fixer.removeHeterogens(keepWater=False)
#
# chains_Todelete = ['a', 'b', 'c', 'd', 'e', 'f', 'A', 'U', 'V', 'W', 'X', 'Y', 'K']
#
# # ALL = ['a', 'b', 'c', 'd', 'e', 'f', 'A', 'U', 'V', 'W', 'X', 'Y', 'K', 'L', 'R']
#
# modeller = Modeller(fixer.topology, fixer.positions)
# # # for r in modeller.topology.chains():
# # #     for i in r.residues():
# # #         print(i.name)
# #
#
#
# toDelete = [r for r in modeller.topology.chains() if r.id in chains_Todelete]
# print(toDelete)
# modeller.delete(toDelete)
# PDBFile.writeFile(
#     modeller.topology,
#     modeller.positions,
#     open(os.path.join('/home/enaz/Desktop/MDPERTOOL_v01/Output', "%s_fixed_ph%s.pdb" % ('pdbid', 7)),
#          "w"),
#     keepIds=True)


# chains = [r.id for r in modeller.topology.chains()]
# app = QtWidgets.QApplication(sys.argv)
# form = ChecklistDialog('Select for Chains Delete', chains, checked=False)
# if form.exec_() == QtWidgets.QDialog.Accepted:
#     print('\n'.join([str(s) for s in form.choices]))
