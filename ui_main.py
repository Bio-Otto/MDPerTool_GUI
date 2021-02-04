import os
import sys
import platform
from PyQt5 import QtCore, QtGui, QtWidgets
#from PyQt5.QtCore import (QPropertyAnimation, QDate, QDateTime, QMetaObject, QObject, QPoint, QRect,
#                          QSize, QTime, QUrl, Qt, QEvent, QRectF)
#from PyQt5.QtGui import (QBrush, QColor, QConicalGradient, QCursor, QFont, QFontDatabase, QIcon, QKeySequence,
#                         QLinearGradient, QPalette, QPainter, QPixmap, QRadialGradient, QRegion)
from PyQt5.QtWidgets import *
import os
import sys
from platform import system, release
from PyQt5 import QtCore, QtGui, QtWidgets
#from PyQt5.QtCore import (QCoreApplication, QPropertyAnimation, QDate, QDateTime, QMetaObject, QObject, QPoint, QRect,
#                          QSize, QTime, QUrl, Qt, QEvent)
# from PyQt5.QtGui import (QBrush, QColor, QConicalGradient, QCursor, QFont, QFontDatabase, QIcon, QKeySequence,
#                          QLinearGradient, QPalette, QPainter, QPixmap, QPainterPath)
# from PyQt5.QtGui import QPainterPath, QRegion
## ==> MAIN WINDOW
# import app_modules
from PyQt5 import uic
from icons import *
import pyqtgraph as pg
# IMPORT FUNCTIONS
from omm_runner import *
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
class MainWindow(QtWidgets.QMainWindow):

    def __init__(self, parent=None):
        super(MainWindow, self).__init__(parent=parent)
        uic.loadUi('MAIN_GUI.ui', self)

        ################################################################################################################
        #                                       ==> START OF WINDOW ATTRIBUTES <==                                     #
        ################################################################################################################
        print('System: ' + platform.system())
        print('Version: ' + platform.release())

        # ------------------------------------ > START OF SIMULATION MONITORING < ------------------------------------ #
        self.created_script = None
        self.Real_Time_Graphs = Graphs()
        self.verticalLayout_16.addWidget(self.Real_Time_Graphs.win)
        # self.setLayout(self.verticalLayout_16)
        # ------------------------------------- > END OF SIMULATION MONITORING < ------------------------------------- #

        # ----- > Remove Standart Title Bar
        UIFunctions.removeTitleBar(True)
        # self.setWindowFlags(QtCore.Qt.FramelessWindowHint | QtCore.Qt.WindowStaysOnTopHint)

        # ----- > Set Window Title
        self.setWindowTitle('MDPerTool v0.1')
        # self.setWindowIcon(QIcon('%s/icons/MDPerTool.ico' % os.getcwd()))
        UIFunctions.labelTitle(self, 'MDPerTool - %s %s' % (platform.system(), platform.release()))
        UIFunctions.labelDescription(self, str(os.getcwd()))

        # ----- > Move The Screen To Center
        qtRectangle = self.frameGeometry()
        centerPoint = QDesktopWidget().availableGeometry().center()
        qtRectangle.moveCenter(centerPoint)
        self.move(qtRectangle.topLeft())

        # ----- > Start With Standart Size
        startSize = QSize(1600, 940)
        self.resize(startSize)
        self.setMinimumSize(startSize)
        # UIFunctions.enableMaximumSize(self, 500, 720)

        # --------------------------------------------- > CREATE MENUS < --------------------------------------------- #
        # ----- > Toggle Menu Size
        self.btn_toggle_menu.clicked.connect(lambda: UIFunctions.toggleMenu(self, 220, True))

        # ----- > Add Custom Menus
        self.stackedWidget.setMinimumWidth(20)
        UIFunctions.addNewMenu(self, "Perturbation", "btn_perturbation", "url(:/20x20/icons/20x20/chemical_20x20.png)",
                               True)
        UIFunctions.addNewMenu(self, "Settings", "btn_settings", "url(:/20x20/icons/20x20/cil-equalizer.png)", False)
        UIFunctions.addNewMenu(self, "About & Contact", "btn_about", "url(:/20x20/icons/20x20/cil-tag.png)", False)

        # ----- > Start Menu Selection
        UIFunctions.selectStandardMenu(self, "btn_perturbation")

        # ----- > Start Page
        self.stackedWidget.setCurrentWidget(self.page_home)

        # ----- > Show / Hide User Icon
        UIFunctions.userIcon(self, "HIO", "", True)

        # ----- > Move Window / Maximize / Restore
        def moveWindow(event):
            # ----- > If Maximized Change to Normal
            if UIFunctions.returStatus() == 1:
                UIFunctions.maximize_restore(self)
            # ----- > Move Window
            if event.buttons() == Qt.LeftButton:
                self.move(self.pos() + event.globalPos() - self.dragPos)
                self.dragPos = event.globalPos()
                event.accept()

        # ----- > Widget to Move
        self.frame_label_top_btns.mouseMoveEvent = moveWindow

        # ----- > Load Definitions
        UIFunctions.uiDefinitions(self)

        ################################################################################################################
        #                                        ==> END OF WINDOW ATTRIBUTES <==                                      #
        ################################################################################################################

        # ------------------------------- > MAINWINDOW WIDGETS FUNCTIONS/PARAMETERS < -------------------------------- #
        self.quit_pushButton.clicked.connect(lambda: UIFunctions.close_application(self))
        self.stop_pushButton.clicked.connect(self.stop_button_clicked)
        self.upload_pdb_Button.clicked.connect(lambda: self.upload_pdb_from_local())
        self.Browse_Output_button.clicked.connect(lambda: self.output_folder_browse())
        self.PDB_ID_lineEdit.textChanged.connect(lambda: Functions.PDB_ID_lineEdit(self))
        self.fetch_pdb_Button.clicked.connect(lambda: self.fetch_and_load_pdbfile())
        self.integrator_kind_comboBox.currentTextChanged.connect(self.Stocasthic_Changed)
        Functions.Send_Available_Platforms_to_GUI(self)
        self.platform_comboBox.currentTextChanged.connect(lambda: Functions.platform_comboBox_Changed(self))
        self.minimize_checkBox.stateChanged.connect(lambda: Functions.minimize_Step_isVisible(self))
        self.State_Data_Reporter_checkBox.stateChanged.connect(lambda: Functions.State_Data_Reporter_Changed(self))
        self.DCD_Reporter_checkBox.stateChanged.connect(lambda: Functions.DCD_Reporter_Changed(self))
        self.XTC_Reporter_checkBox.stateChanged.connect(lambda: Functions.XTC_Reporter_Changed(self))
        self.Run.clicked.connect(self.run_btn_clicked)

        # ----------------------------------- > START OF PYMOL RELEATED BUTTONS < ------------------------------------ #
        UIFunctions.start_pymol(self)
        self.add_residue_pushButton.clicked.connect(lambda: Functions.add_residue_toList(self))
        self.discard_residue_pushButton.clicked.connect(lambda: Functions.discard_residue_fromList(self))
        self.selected_residues_listWidget.itemDoubleClicked.connect(lambda: UIFunctions.show_residue_labels(self))
        self.refresh_pushButton.clicked.connect(lambda: UIFunctions.clear_residue_labels(self))
        self.activate_pymol_navigation.clicked.connect(lambda: UIFunctions.activate_navigation_on_Pymol(self))
        self.deactivate_pymol_navigation.clicked.connect(lambda: UIFunctions.deactivate_navigation_on_Pymol(self))
        self.visualization_Handel_buttons_changing()
        self.Handel_Buttons()
        self.ss_beatiful_snapshoot.clicked.connect(lambda: UIFunctions.show_beatiful_in_Pymol(self))
        self.get_figure_pushButton.clicked.connect(lambda: UIFunctions.save_as_png_Pymol(self))
        self.Handel_Save_Figure_Options_Changed()
        self.Handel_Save_Figure_Options()

    def run_btn_clicked(self):
        self.Run.setEnabled(False)
        self.start_monitoring = False
        self.start_monitoring = Advanced.send_arg_to_Engine(self)

        if self.start_monitoring:
            self.show_simulation_monitoring()

        if not self.start_monitoring:
            self.Run.setEnabled(False)

    def stop_button_clicked(self):
        self.__stop = True
        try:
            self.Real_Time_Graphs.stop_th()

        except Exception as ins:
            QMessageBox.warning(self, "The program can't stop the running Simulation", str(ins))

    def show_simulation_monitoring(self):
        self.stackedWidget.setCurrentIndex(1)
        self.Real_Time_Graphs.run_script(self.created_script)

    def fetch_and_load_pdbfile(self):
        """
        Fetch action for getting crystal structure from pdb databank
        :return: path of downloaded pdb file or None conditions will return
        """
        fetched_and_modified_pdb = Functions.Fetch_PDB_File(self)

        if fetched_and_modified_pdb:
            UIFunctions.load_pdb_to_pymol(self, fetched_and_modified_pdb)

        if fetched_and_modified_pdb is None:
            Message_Boxes.Information_message(self, 'Wrong pdb id "%s"' % self.PDB_ID_lineEdit.text(),
                                              'There is no protein crystal structure in this expression.',
                                              Style.MessageBox_stylesheet)

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
                pdb_fix_dialog_answer = checked_list.exec_()
                if pdb_fix_dialog_answer == QtWidgets.QDialog.Accepted:
                    selected_chains = [str(s) for s in checked_list.choices]

                    delete_chains = list(set(chains) - set(selected_chains))

                    modified_pdb = pdb_Tools.fetched_pdb_fix(self, pdb_path, self.Output_Folder_textEdit.toPlainText(),
                                                             ph=7, chains_to_remove=delete_chains)

                    self.upload_pdb_textEdit.setText(modified_pdb)

                    self.combobox = Helper_Functions.fill_residue_combobox(self, modified_pdb)
                    for i in self.combobox:
                        self.res1_comboBox.addItem(str(i))
                    self.res1_comboBox.clear()  # delete all items from comboBox
                    self.res1_comboBox.addItems(self.combobox)  # add the actual content of self.comboData

                elif pdb_fix_dialog_answer == QtWidgets.QDialog.Rejected:
                    modified_pdb = pdb_Tools.fetched_pdb_fix(self, pdb_path, self.Output_Folder_textEdit.toPlainText(),
                                                             ph=7, chains_to_remove=None)

                    self.upload_pdb_textEdit.setText(modified_pdb)

                    self.combobox = Helper_Functions.fill_residue_combobox(self, modified_pdb)
                    for i in self.combobox:
                        self.res1_comboBox.addItem(str(i))
                    self.res1_comboBox.clear()  # delete all items from comboBox
                    self.res1_comboBox.addItems(self.combobox)  # add the actual content of self.comboData

                UIFunctions.load_pdb_to_pymol(self, modified_pdb)

        except TypeError:
            pass

    def output_folder_browse(self):
        Functions.output_file(self)

    def Stocasthic_Changed(self):
        Functions.Stochastic_changed(self)

    ####################################################################################################################
    #                                     ==> START OF DYNAMIC MENUS FUNCTIONS < ==                                    #
    ####################################################################################################################
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

    ####################################################################################################################
    #                                     == > END OF DYNAMIC MENUS FUNCTIONS < ==                                     #
    ####################################################################################################################
    ####################################################################################################################
    #                                          == > START OF APP EVENTS < ==                                           #
    ####################################################################################################################

    # ----- > Mouse Double Click
    def eventFilter(self, watched, event):
        if watched == self.le and event.type() == QtCore.QEvent.MouseButtonDblClick:
            print("pos: ", event.pos())

    # ----- > Mouse Click Event
    def mousePressEvent(self, event):
        self.dragPos = event.globalPos()
        if event.buttons() == Qt.LeftButton:
            print('Mouse click: LEFT CLICK')
        if event.buttons() == Qt.RightButton:
            print('Mouse click: RIGHT CLICK')
        if event.buttons() == Qt.MidButton:
            print('Mouse click: MIDDLE BUTTON')

    # ----- > Key Pressed
    def keyPressEvent(self, event):
        print('Key: ' + str(event.key()) + ' | Text Press: ' + str(event.text()))

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
        print('Height: ' + str(self.height()) + ' | Width: ' + str(self.width()))

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
    ####################################################################################################################
    #                                       == > START OF SPLASH SCREEN < ==                                           #
    ####################################################################################################################


class SplashScreen(QMainWindow):

    def __init__(self, parent=None, *args, **kwargs):
        super().__init__(parent=parent, *args, **kwargs)

        uic.loadUi('splash_screen.ui', self)

        # ----- > Set App Icon
        app_icon = QtGui.QIcon()
        app_icon.addFile('%s/icons/big_icons/style_icon_48x48.png' % os.getcwd())

        # ----- > Move to The Screen Center
        qtRectangle = self.frameGeometry()
        centerPoint = QDesktopWidget().availableGeometry().center()
        qtRectangle.moveCenter(centerPoint)
        self.move(qtRectangle.topLeft())

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
        self.timer.start(15)

        # ----- > Change Description on Splash Screen
        self.label_title.setText("MDPerTool v0.1")
        self.label_credits.setText("<strong>Devoloped</strong> by Ozbek's Lab")

        # ----- > Initial Text on Splash Screen
        self.label_description.setText("<strong>WELCOME</strong> TO MDPerTool V0.1 PLATFORM")

        # ----- > CHANGE TEXT ON SPLASH SCREEN
        QtCore.QTimer.singleShot(300, lambda: self.label_description.setText("<strong>LOADING</strong> ENVIRONMENT"))
        QtCore.QTimer.singleShot(250, lambda: self.label_description.setText("<strong>LOADING</strong> USER INTERFACE"))

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
            self.main.show()
            self.close()

        counter += 1
    # ------------------------------------------- > END OF APP FUNCTIONS < ------------------------------------------- #


if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    app.setStyleSheet(Style.QToolTip_stylesheet)
    window = SplashScreen()
    sys.exit(app.exec_())
