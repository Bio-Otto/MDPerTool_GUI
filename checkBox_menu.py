
from PySide2 import QtCore, QtWidgets, QtGui
from PySide2.QtGui import QIcon,QFont,QPixmap,QPalette
from PySide2.QtCore import QCoreApplication, Qt,QBasicTimer, QPoint
from PySide2.QtWidgets import QDesktopWidget

stylesheet = """

        QWidget
        {
            color: #b1b1b1;
            background-color: #323232;
            selection-background-color:#323232;
            selection-color: black;
            background-clip: border;
            border-image: none;
            border: 0px transparent black;
            outline: 0;
            border: 1.5px solid;
            border-radius: 5px;
            padding: 8px;
        }
        
        QLineEdit 
        {
            background-color: rgb(27, 29, 35); 
            border-radius: 5px;
        }
        
        QLineEdit:hover 
        {
            border: 2px solid rgb(64, 71, 88);
        }
        
        QLineEdit:focus 
        {
            border: 2px solid rgb(91, 101, 124);
        }
        
        QCheckBox 
        {
            spacing: 5px;
        } 
        
        QCheckBox::indicator 
        {
            width: 15px; 
            height: 15px;
        }
        
        QCheckBox::indicator:unchecked 
        {
            border-image: url(:/16x16/icons/16x16/cil-circle.png);
        }
        
        QCheckBox::indicator:checked 
        {
            border-image: url(:/20x20/icons/20x20/cil-check.png);
        }
        
        QListView 
        {
            font-size: 12pt; 
            background : rgb(64, 71, 88); 
            border: 2px solid grey; 
            border-radius: 5px; 
            text-align: center; 
            font-weight: bold;
            margin-top:10px;
            margin-left:3px;
            margin-right:2px;
        }
        
        QListView::item:!selected:hover 
        {
            background: rgb(64, 71, 88); 
            outline: 0; 
            color: #eff0f1; 
            background-color: rgb(15, 133, 163);
        }
        
        QListView::item:selected:hover
        {
        background : rgb(64, 71, 88);
        }
                
        QPushButton 
        {
            color: white; 
            font-weight: bold; 
            font-size: 10px; 
            border: 2px solid rgb(52, 59, 72); 
            border-radius: 5px; 
            background-color:  rgb(22, 200, 244); 
            margin-top:1px; 
            margin-bottom: 1px; 
            border - width: 1px; 
            padding: 5px; 
            outline: none;
        }
       
        QPushButton:hover 
        { 
            background-color: rgb(255, 17, 100); 
            border: 2px solid rgb(61, 70, 86);
        }
        
        QPushButton:pressed 
        { 
            background-color:  rgb(15, 133, 163); 
            border: 2px solid rgb(43, 50, 61);
        }

        QScrollBar:vertical
        {
            background-color: #323232;
            width: 15px;
            margin: 15px 3px 15px 3px;
            border: 1px transparent #2A2929;
        }
        
        QScrollBar::handle:vertical
        {
            background-color: rgb(255, 17, 100);
            min-height: 5px;
            border-radius: 3px;
            
        }
        
        QScrollBar::sub-line:vertical
        {
            border: none;
            background: none;
            color: none;
        }
        
        QScrollBar::add-line:vertical
        {
            border: none;
            background: none;
            color: none;
        }
        
        QScrollBar::add-page:vertical, QScrollBar::sub-page:vertical
        {
            background: none;
        }

    """


class ChecklistDialog(QtWidgets.QDialog):

    def __init__(self,
                 name,
                 stringlist=None,
                 checked=False,
                 icon=None,
                 parent=None,
                 ):
        super(ChecklistDialog, self).__init__(parent)
        self.setWindowFlags(QtCore.Qt.Dialog | QtCore.Qt.FramelessWindowHint)
        self.setWindowFlags(QtCore.Qt.FramelessWindowHint | QtCore.Qt.WindowStaysOnTopHint)
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
        self.cancelButton = QtWidgets.QPushButton('Don\'t Fix')
        self.selectButton = QtWidgets.QPushButton('Select All')
        self.unselectButton = QtWidgets.QPushButton('Unselect All')

        hbox = QtWidgets.QHBoxLayout()
        # hbox.addStretch(1)
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

        # colors = ['#7fc97f', '#beaed4', '#fdc086', '#ffff99', '#386cb0', '#f0027f', '#bf5b17', '#666666']
        # self.listView.setStyleSheet("QListWidget::item { border-bottom: 1px solid red; background-color: blue;}")

        self.setStyleSheet(stylesheet)
        # center
        self.oldPos = self.pos()

    def center(self):
        qr = self.frameGeometry()
        cp = QDesktopWidget().availableGeometry().center()
        qr.moveCenter(cp)
        self.move(qr.topLeft())

    def mousePressEvent(self, event):
        self.oldPos = event.globalPos()

    def mouseMoveEvent(self, event):
        delta = QPoint(event.globalPos() - self.oldPos)
        # print(delta)
        self.move(self.x() + delta.x(), self.y() + delta.y())
        self.oldPos = event.globalPos()

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
#     'Banana',
#     'Apple',
#     'Elderberry',
#     'Clementine',
#     'Fig',
#     'Guava',
#     'Mango',
#     'Honeydew Melon'
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
#
# form = ChecklistDialog('Select for Chains Delete', fruits, checked=False)
# if form.exec_() == QtWidgets.QDialog.Accepted:
#     print('\n'.join([str(s) for s in form.choices]))
