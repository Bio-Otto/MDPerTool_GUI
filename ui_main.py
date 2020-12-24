import os
import sys
import platform
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtCore import (QCoreApplication, QPropertyAnimation, QDate, QDateTime, QMetaObject, QObject, QPoint, QRect,
                          QSize, QTime, QUrl, Qt, QEvent, QRectF)
from PyQt5.QtGui import (QBrush, QColor, QConicalGradient, QCursor, QFont, QFontDatabase, QIcon, QKeySequence,
                         QLinearGradient, QPalette, QPainter, QPixmap, QRadialGradient, QRegion)
from PyQt5.QtWidgets import *
import os
import sys
from platform import system, release
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtCore import (QCoreApplication, QPropertyAnimation, QDate, QDateTime, QMetaObject, QObject, QPoint, QRect,
                          QSize, QTime, QUrl, Qt, QEvent)
from PyQt5.QtGui import (QBrush, QColor, QConicalGradient, QCursor, QFont, QFontDatabase, QIcon, QKeySequence,
                         QLinearGradient, QPalette, QPainter, QPixmap, QPainterPath)

## ==> MAIN WINDOW
# import app_modules
from PyQt5 import uic
from icons import *
import pyqtgraph as pg
# IMPORT FUNCTIONS
# from omm_runner import *
from ui_functions import *
from builder import *
from omm_runner import *
from pdbfixer import PDBFixer

# from app_functions import *
# from browse_menu import browse_file
# from omm_runner import *
## ==> GLOBALS
counter = 0



# YOUR APPLICATION
class MainWindow(QMainWindow):
    def __init__(self, parent=None):
        super(MainWindow, self).__init__(parent=parent)
        uic.loadUi('MAIN_GUI.ui', self)

        ## PRINT ==> SYSTEM
        print('System: ' + platform.system())
        print('Version: ' + platform.release())

        ########################################################################
        ## START - WINDOW ATTRIBUTES
        ########################################################################

        ## START - SIMULATION MONITORING
        self.created_script = None
        self.Real_Time_Graphs = Graphs()
        self.verticalLayout_16.addWidget(self.Real_Time_Graphs.win)
        self.setLayout(self.verticalLayout_16)
        ## END - SIMULATION MONITORING

        ## REMOVE ==> STANDARD TITLE BAR

        UIFunctions.start_pymol(self)
        UIFunctions.removeTitleBar(True)

        # self.setWindowFlags(QtCore.Qt.FramelessWindowHint | QtCore.Qt.WindowStaysOnTopHint)


        ## ==> END ##

        ## SET ==> WINDOW TITLE
        self.setWindowTitle('MDPerTool v0.1')
        # self.setWindowIcon(QIcon('%s/icons/MDPerTool.ico' % os.getcwd()))
        UIFunctions.labelTitle(self, 'MDPerTool - %s %s' % (platform.system(), platform.release()))
        UIFunctions.labelDescription(self, str(os.getcwd()))
        ## ==> END ##

        ## MOVE TO THE SCREEN CENTER
        qtRectangle = self.frameGeometry()
        centerPoint = QDesktopWidget().availableGeometry().center()
        qtRectangle.moveCenter(centerPoint)
        self.move(qtRectangle.topLeft())
        ## ==> END ##

        ## REMOVE ==> STANDARD TITLE BAR
        startSize = QSize(1600, 940)
        self.resize(startSize)
        self.setMinimumSize(startSize)
        # UIFunctions.enableMaximumSize(self, 500, 720)
        ## ==> END ##

        ## ==> CREATE MENUS
        ########################################################################

        ## ==> TOGGLE MENU SIZE
        self.btn_toggle_menu.clicked.connect(lambda: UIFunctions.toggleMenu(self, 220, True))
        ## ==> END ##

        ## ==> ADD CUSTOM MENUS
        self.stackedWidget.setMinimumWidth(20)
        UIFunctions.addNewMenu(self, "Perturbation", "btn_perturbation", "url(:/20x20/icons/20x20/chemical_20x20.png)",
                               True)
        UIFunctions.addNewMenu(self, "Settings", "btn_settings", "url(:/20x20/icons/20x20/cil-equalizer.png)", False)
        UIFunctions.addNewMenu(self, "About & Contact", "btn_about", "url(:/20x20/icons/20x20/cil-tag.png)", False)
        ## ==> END ##

        # START MENU => SELECTION
        UIFunctions.selectStandardMenu(self, "btn_home")
        ## ==> END ##

        ## ==> START PAGE
        self.stackedWidget.setCurrentWidget(self.page_home)
        ## ==> END ##

        ## USER ICON ==> SHOW HIDE
        UIFunctions.userIcon(self, "HIO", "", True)

        ## ==> END ##

        ## ==> MOVE WINDOW / MAXIMIZE / RESTORE
        ########################################################################
        def moveWindow(event):
            # IF MAXIMIZED CHANGE TO NORMAL
            if UIFunctions.returStatus() == 1:
                UIFunctions.maximize_restore(self)

            # MOVE WINDOW
            if event.buttons() == Qt.LeftButton:
                self.move(self.pos() + event.globalPos() - self.dragPos)
                self.dragPos = event.globalPos()
                event.accept()

        # WIDGET TO MOVE
        self.frame_label_top_btns.mouseMoveEvent = moveWindow
        ## ==> END ##

        ## ==> LOAD DEFINITIONS
        ########################################################################
        UIFunctions.uiDefinitions(self)
        ## ==> END ##

        ########################################################################
        ## END - WINDOW ATTRIBUTES
        ############################## ---/--/--- ##############################

        ########################################################################
        #                                                                      #
        ## START -------------- WIDGETS FUNCTIONS/PARAMETERS ---------------- ##
        #                                                                      #
        ## ==> USER CODES BELLOW                                              ##
        ########################################################################

        self.upload_pdb_Button.clicked.connect(lambda: self.upload_pdb_from_local())
        self.Browse_Output_button.clicked.connect(lambda: self.output_folder_browse())
        self.PDB_ID_lineEdit.textChanged.connect(lambda: Functions.PDB_ID_lineEdit(self))
        self.fetch_pdb_Button.clicked.connect(lambda: Functions.Fetch_PDB_File(self))

        self.integrator_kind_comboBox.currentTextChanged.connect(self.Stocasthic_Changed)
        Functions.Send_Available_Platforms_to_GUI(self)
        self.platform_comboBox.currentTextChanged.connect(lambda: Functions.platform_comboBox_Changed(self))
        self.minimize_checkBox.stateChanged.connect(lambda: Functions.minimize_Step_isVisible(self))
        self.State_Data_Reporter_checkBox.stateChanged.connect(lambda: Functions.State_Data_Reporter_Changed(self))
        self.DCD_Reporter_checkBox.stateChanged.connect(lambda: Functions.DCD_Reporter_Changed(self))
        self.XTC_Reporter_checkBox.stateChanged.connect(lambda: Functions.XTC_Reporter_Changed(self))

        self.Run.clicked.connect(self.run_btn_clicked)

        ########################################################################
        #
        ## END --------------- WIDGETS FUNCTIONS/PARAMETERS ----------------- ##
        #                                                                      #
        ############################## ---/--/--- ##############################

        ## SHOW ==> MAIN WINDOW
        ########################################################################

        ## ==> END ##
    def run_btn_clicked(self):
        self.start_monitoring = Advanced.send_arg_to_Engine(self)
        if self.start_monitoring:
            self.show_simulation_monitoring()

    def show_simulation_monitoring(self):
        self.stackedWidget.setCurrentIndex(1)
        self.Real_Time_Graphs.run_script(self.created_script)


    def upload_pdb_from_local(self):
        global selected_chains
        try:
            upload_condition, pdb_path = Functions.browse_pdbFile(self)

            if os.path.exists(pdb_path):
                fixer = PDBFixer(pdb_path)
                fixer.removeHeterogens(keepWater=False)

                modeller = Modeller(fixer.topology, fixer.positions)
                chains = [r.id for r in modeller.topology.chains()]

                checked_list = ChecklistDialog('Select the chain (s) to be used in the system', chains,
                                               checked=True)
                if checked_list.exec_() == QtWidgets.QDialog.Accepted:
                    selected_chains = [str(s) for s in checked_list.choices]

                delete_chains = list(set(chains) - set(selected_chains))

                modified_pdb = pdb_Tools.fetched_pdb_fix(self, pdb_path, self.Output_Folder_textEdit.toPlainText(),
                                                         ph=7, chains_to_remove=delete_chains)

                self.upload_pdb_textEdit.setText(modified_pdb)

                self.combobox = Helper_Functions.fill_residue_combobox(self, modified_pdb)
                for i in self.combobox:
                    self.res1_comboBox.addItem(str(i))
                    self.res2_comboBox.addItem(str(i))
                self.res1_comboBox.clear()  # delete all items from comboBox
                self.res1_comboBox.addItems(self.combobox)  # add the actual content of self.comboData
                self.res2_comboBox.clear()  # delete all items from comboBox
                self.res2_comboBox.addItems(self.combobox)  # add the actual content of self.comboData

        except Exception as instance:
            pass
            # QMessageBox.critical(self, 'An error occurred while loading files.', repr(instance))

    def output_folder_browse(self):
        Functions.output_file(self)

    def Stocasthic_Changed(self):
        Functions.Stochastic_changed(self)

    ########################################################################
    ## MENUS ==> DYNAMIC MENUS FUNCTIONS
    ########################################################################
    def Button(self):
        # GET BT CLICKED
        btnWidget = self.sender()

        # PAGE HOME
        if btnWidget.objectName() == "btn_home":
            self.stackedWidget.setCurrentWidget(self.page_home)
            UIFunctions.resetStyle(self, "btn_home")
            UIFunctions.labelPage(self, "Home")
            btnWidget.setStyleSheet(UIFunctions.selectMenu(btnWidget.styleSheet()))

        # PAGE PERTURBATION
        if btnWidget.objectName() == "btn_perturbation":
            self.stackedWidget.setCurrentWidget(self.page_home)
            UIFunctions.resetStyle(self, "btn_perturbation")
            UIFunctions.labelPage(self, "Perturbation User Interface")
            btnWidget.setStyleSheet(UIFunctions.selectMenu(btnWidget.styleSheet()))

        # PAGE WIDGETS
        if btnWidget.objectName() == "btn_settings":
            self.stackedWidget.setCurrentWidget(self.page_settings)
            UIFunctions.resetStyle(self, "btn_settings")
            UIFunctions.labelPage(self, "Settings")
            btnWidget.setStyleSheet(UIFunctions.selectMenu(btnWidget.styleSheet()))

        # PAGE ABOUT
        if btnWidget.objectName() == "btn_about":
            self.stackedWidget.setCurrentWidget(self.page_contact)
            UIFunctions.resetStyle(self, "btn_about")
            UIFunctions.labelPage(self, "ABOUT & CONTACT")
            btnWidget.setStyleSheet(UIFunctions.selectMenu(btnWidget.styleSheet()))

    ## ==> END ##

    ########################################################################
    ## START ==> APP EVENTS
    ########################################################################

    ## EVENT ==> MOUSE DOUBLE CLICK
    ########################################################################
    def eventFilter(self, watched, event):
        if watched == self.le and event.type() == QtCore.QEvent.MouseButtonDblClick:
            print("pos: ", event.pos())

    ## ==> END ##

    ## EVENT ==> MOUSE CLICK
    ########################################################################
    def mousePressEvent(self, event):
        self.dragPos = event.globalPos()
        if event.buttons() == Qt.LeftButton:
            print('Mouse click: LEFT CLICK')
        if event.buttons() == Qt.RightButton:
            print('Mouse click: RIGHT CLICK')
        if event.buttons() == Qt.MidButton:
            print('Mouse click: MIDDLE BUTTON')

    ## ==> END ##

    ## EVENT ==> KEY PRESSED
    ########################################################################
    def keyPressEvent(self, event):
        print('Key: ' + str(event.key()) + ' | Text Press: ' + str(event.text()))

    ## ==> END ##

    ## EVENT ==> RESIZE EVENT
    ########################################################################
    def resizeEvent(self, event):
        self.resizeFunction()
        path = QPainterPath()
        path.addRoundedRect(QRectF(self.rect()), 5, 5)
        reg = QRegion(path.toFillPolygon().toPolygon())
        self.setMask(reg)
        return super(MainWindow, self).resizeEvent(event)

    def resizeFunction(self):
        print('Height: ' + str(self.height()) + ' | Width: ' + str(self.width()))
    ## ==> END ##



# SPLASH SCREEN
class SplashScreen(QMainWindow):

    def __init__(self):
        super().__init__()
        uic.loadUi('splash_screen.ui', self)

        ## UI ==> INTERFACE CODES
        # set app icon
        app_icon = QtGui.QIcon()
        app_icon.addFile('%s/icons/big_icons/style_icon_48x48.png' % os.getcwd())
        # app_icon.addFile('gui/icons/24x24.png', QtCore.QSize(24, 24))
        # app_icon.addFile('gui/icons/32x32.png', QtCore.QSize(32, 32))
        # app_icon.addFile('gui/icons/48x48.png', QtCore.QSize(48, 48))
        # app_icon.addFile('gui/icons/256x256.png', QtCore.QSize(256, 256))

        ## MOVE TO THE SCREEN CENTER
        qtRectangle = self.frameGeometry()
        centerPoint = QDesktopWidget().availableGeometry().center()
        qtRectangle.moveCenter(centerPoint)
        self.move(qtRectangle.topLeft())
        ## ==> END ##

        ########################################################################
        self.setWindowIcon(QIcon(app_icon))
        ## REMOVE TITLE BAR
        self.setWindowFlags(QtCore.Qt.FramelessWindowHint | QtCore.Qt.WindowStaysOnTopHint)
        self.setAttribute(QtCore.Qt.WA_TranslucentBackground)

        ## DROP SHADOW EFFECT
        self.shadow = QGraphicsDropShadowEffect(self)
        self.shadow.setBlurRadius(20)
        self.shadow.setXOffset(0)
        self.shadow.setYOffset(0)
        self.shadow.setColor(QColor(0, 0, 0, 60))
        self.dropShadowFrame.setGraphicsEffect(self.shadow)

        ## QTIMER ==> START
        self.timer = QtCore.QTimer()
        self.timer.timeout.connect(self.progress)
        # TIMER IN MILLISECONDS
        self.timer.start(5)

        # CHANGE DESCRIPTION
        self.label_title.setText("MDPerTool v0.1")
        self.label_credits.setText("<strong>Devoloped</strong> by Ozbek's Lab")
        # Initial Text
        self.label_description.setText("<strong>WELCOME</strong> TO MDPerTool V0.1 PLATFORM")

        # Change Texts
        QtCore.QTimer.singleShot(100, lambda: self.label_description.setText("<strong>LOADING</strong> ENVIRONMENT"))
        QtCore.QTimer.singleShot(200,
                                 lambda: self.label_description.setText("<strong>LOADING</strong> USER INTERFACE"))

        ## SHOW ==> MAIN WINDOW
        ########################################################################
        self.show()
        ## ==> END ##

    ## ==> APP FUNCTIONS
    ########################################################################
    def progress(self):
        global counter

        # SET VALUE TO PROGRESS BAR
        self.progressBar.setValue(counter)

        # CLOSE SPLASH SCREE AND OPEN APP
        if counter > 100:
            # STOP TIMER
            self.timer.stop()

            # SHOW MAIN WINDOW
            self.main = MainWindow()
            self.main.show()

            # CLOSE SPLASH SCREEN
            self.close()

        # INCREASE COUNTER
        counter += 1


if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = SplashScreen()
    sys.exit(app.exec_())
