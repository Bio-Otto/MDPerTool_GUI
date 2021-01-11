from PyQt5.QtWidgets import *
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtCore import (QCoreApplication, QPropertyAnimation, QDate, QDateTime, QMetaObject, QObject, QPoint, QRect,
                          QSize, QTime, QUrl, Qt, QEvent)
from PyQt5.QtGui import (QBrush, QColor, QConicalGradient, QCursor, QFont, QFontDatabase, QIcon, QKeySequence,
                         QLinearGradient, QPalette, QPainter, QPixmap, QRadialGradient)

## ==> GUI FILE
from ui_main import *
# IMPORT QSS CUSTOM
from ui_styles import Style
import sip

## ==> GLOBALS
GLOBAL_STATE = 0
GLOBAL_TITLE_BAR = True

## ==> COUT INITIAL MENU
count = 1
from app_functions import *
from PyMolWidget import PymolQtWidget


class UIFunctions(MainWindow):
    ## ==> GLOBALS
    GLOBAL_STATE = 0
    GLOBAL_TITLE_BAR = True

    ########################################################################
    ## START - GUI FUNCTIONS
    ########################################################################

    ## ==> CLOSE APPLICATION
    @staticmethod
    def close_application(self):
        try:
            """Close Application Question Message Box."""
            close_program_msgbox = QMessageBox(QMessageBox.Question, "Be carefull !",
                                               "Do you really want to close program?")

            close_program_msgbox.setIcon(QMessageBox.Question)
            close_program_msgbox.addButton(QMessageBox.Yes)
            close_program_msgbox.addButton(QMessageBox.No)
            close_program_msgbox.setDefaultButton(QMessageBox.No)
            # msgbox.setCheckBox(cb)
            # msg.setWindowIcon(QIcon(ICON_PATH))
            close_program_msgbox.setStyleSheet(Style.MessageBox_stylesheet)
            close_program_msgbox.setWindowFlags(QtCore.Qt.FramelessWindowHint | QtCore.Qt.WindowStaysOnTopHint)
            close_answer = close_program_msgbox.exec_()

            if close_answer == QMessageBox.Yes:
                self.close()

            if close_answer == QMessageBox.No:
                pass

        except Exception as inst:
            print(inst)

    ## ==> MAXIMIZE/RESTORE
    ########################################################################
    def maximize_restore(self):
        global GLOBAL_STATE
        status = GLOBAL_STATE
        if status == 0:
            self.showMaximized()
            GLOBAL_STATE = 1
            self.horizontalLayout.setContentsMargins(0, 0, 0, 0)
            self.btn_maximize_restore.setToolTip("Restore")
            self.btn_maximize_restore.setIcon(QtGui.QIcon(u":/16x16/icons/16x16/cil-window-restore.png"))
            self.frame_top_btns.setStyleSheet("background-color: rgb(27, 29, 35)")
            self.frame_size_grip.hide()
        else:
            GLOBAL_STATE = 0
            self.showNormal()
            self.resize(self.width() + 1, self.height() + 1)
            self.horizontalLayout.setContentsMargins(10, 10, 10, 10)
            self.btn_maximize_restore.setToolTip("Maximize")
            self.btn_maximize_restore.setIcon(QtGui.QIcon(u":/16x16/icons/16x16/cil-window-maximize.png"))
            self.frame_top_btns.setStyleSheet("background-color: rgba(27, 29, 35, 200)")
            self.frame_size_grip.show()

    ## ==> RETURN STATUS
    @staticmethod
    def returStatus():
        return GLOBAL_STATE

    ## ==> SET STATUS
    @staticmethod
    def setStatus(status):
        global GLOBAL_STATE
        GLOBAL_STATE = status

    ## ==> ENABLE MAXIMUM SIZE
    ########################################################################
    def enableMaximumSize(self, width, height):
        if width != '' and height != '':
            self.setMaximumSize(QSize(width, height))
            self.frame_size_grip.hide()
            self.btn_maximize_restore.hide()

    ## ==> TOGGLE MENU
    ########################################################################
    def toggleMenu(self, maxWidth, enable):
        if enable:
            # GET WIDTH
            width = self.frame_left_menu.width()
            maxExtend = maxWidth
            standard = 70

            # SET MAX WIDTH
            if width == 70:
                widthExtended = maxExtend
            else:
                widthExtended = standard

            # ANIMATION
            self.animation = QPropertyAnimation(self.frame_left_menu, b"minimumWidth")
            self.animation.setDuration(300)
            self.animation.setStartValue(width)
            self.animation.setEndValue(widthExtended)
            self.animation.setEasingCurve(QtCore.QEasingCurve.InOutQuart)
            self.animation.start()

    ## ==> SET TITLE BAR
    ########################################################################
    def removeTitleBar(self):
        global GLOBAL_TITLE_BAR
        GLOBAL_TITLE_BAR = self

    ## ==> HEADER TEXTS
    ########################################################################
    # LABEL TITLE
    def labelTitle(self, text):
        self.label_title_bar_top.setText(text)

    # LABEL DESCRIPTION
    def labelDescription(self, text):
        self.label_top_info_1.setText(text)

    ## ==> DYNAMIC MENUS
    ########################################################################
    def addNewMenu(self, name, objName, icon, isTopMenu):
        font = QFont()
        font.setFamily(u"Segoe UI")
        button = QPushButton(str(count), self)
        button.setObjectName(objName)
        sizePolicy3 = QSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        sizePolicy3.setHorizontalStretch(0)
        sizePolicy3.setVerticalStretch(0)
        sizePolicy3.setHeightForWidth(button.sizePolicy().hasHeightForWidth())
        button.setSizePolicy(sizePolicy3)
        button.setMinimumSize(QSize(0, 70))
        button.setLayoutDirection(Qt.LeftToRight)
        button.setFont(font)
        button.setStyleSheet(Style.style_bt_standard.replace('ICON_REPLACE', icon))
        button.setText(name)
        button.setToolTip(name)
        button.clicked.connect(self.Button)

        if isTopMenu:
            self.layout_menus.addWidget(button)
        else:
            self.layout_menu_bottom.addWidget(button)

    ## ==> SELECT/DESELECT MENU
    ########################################################################
    ## ==> SELECT
    def selectMenu(getStyle):
        select = getStyle + ("QPushButton { border-right: 7px solid rgb(44, 49, 60); }")
        return select

    ## ==> DESELECT
    def deselectMenu(getStyle):
        deselect = getStyle.replace("QPushButton { border-right: 7px solid rgb(44, 49, 60); }", "")
        return deselect

    ## ==> START SELECTION
    def selectStandardMenu(self, widget):
        for w in self.frame_left_menu.findChildren(QPushButton):
            if w.objectName() == widget:
                w.setStyleSheet(UIFunctions.selectMenu(w.styleSheet()))

    ## ==> RESET SELECTION
    def resetStyle(self, widget):
        for w in self.frame_left_menu.findChildren(QPushButton):
            if w.objectName() != widget:
                w.setStyleSheet(UIFunctions.deselectMenu(w.styleSheet()))

    ## ==> CHANGE PAGE LABEL TEXT
    def labelPage(self, text):
        newText = '| ' + text.upper()
        self.label_top_info_2.setText(newText)

    ## ==> USER ICON
    ########################################################################
    def userIcon(self, initialsTooltip, icon, showHide):
        if showHide:
            # SET TEXT
            self.label_user_icon.setText(initialsTooltip)

            # SET ICON
            if icon:
                style = self.label_user_icon.styleSheet()
                setIcon = "QLabel { background-image: " + icon + "; }"
                self.label_user_icon.setStyleSheet(style + setIcon)
                self.label_user_icon.setText('')
                self.label_user_icon.setToolTip(initialsTooltip)
        else:
            self.label_user_icon.hide()

    ########################################################################
    ## END - GUI FUNCTIONS
    ########################################################################

    ########################################################################
    ## START - GUI DEFINITIONS
    ########################################################################

    ## ==> UI DEFINITIONS
    ########################################################################
    def uiDefinitions(self):
        def dobleClickMaximizeRestore(event):
            # IF DOUBLE CLICK CHANGE STATUS
            if event.type() == QtCore.QEvent.MouseButtonDblClick:
                QtCore.QTimer.singleShot(250, lambda: UIFunctions.maximize_restore(self))

        ## REMOVE ==> STANDARD TITLE BAR
        if GLOBAL_TITLE_BAR:
            self.setStyleSheet("background:rgb(27, 29, 35);")
            self.setWindowFlags(QtCore.Qt.FramelessWindowHint)
            self.frame_main.setAttribute(QtCore.Qt.WA_TranslucentBackground)
            # self.centralwidget.setAttribute(QtCore.Qt.WA_NoSystemBackground)
            # self.setStyleSheet("background:rgb(27, 29, 35);")
            self.frame_label_top_btns.mouseDoubleClickEvent = dobleClickMaximizeRestore
        else:
            self.horizontalLayout.setContentsMargins(0, 0, 0, 0)
            self.frame_label_top_btns.setContentsMargins(8, 0, 0, 5)
            self.frame_label_top_btns.setMinimumHeight(42)
            self.frame_icon_top_bar.hide()
            self.frame_btns_right.hide()
            self.frame_size_grip.hide()

        ## SHOW ==> DROP SHADOW
        # self.shadow = QGraphicsDropShadowEffect(self)
        # self.shadow.setBlurRadius(17)
        # self.shadow.setXOffset(0)
        # self.shadow.setYOffset(0)
        # self.shadow.setColor(QColor(0, 0, 0, 50))
        # self.frame_main.setGraphicsEffect(self.shadow)
        ## DROP SHADOW EFFECT

        ## ==> RESIZE WINDOW
        self.sizegrip = QSizeGrip(self.frame_size_grip)
        self.sizegrip.setStyleSheet("width: 20px; height: 20px; margin 0px; padding: 0px;")

        ### ==> MINIMIZE
        self.btn_minimize.clicked.connect(lambda: self.showMinimized())

        ## ==> MAXIMIZE/RESTORE
        self.btn_maximize_restore.clicked.connect(lambda: UIFunctions.maximize_restore(self))

        ## SHOW ==> CLOSE APPLICATION
        self.btn_close.clicked.connect(lambda: self.close())

    ########################################################################
    ## END - GUI DEFINITIONS
    ########################################################################

    ########################################################################
    ## == > OPEN SOURCE PYMOL 2.4 INTEGRATION START
    # def start_pymol(self):
    #     self.ProteinView = PymolQtWidget(self)
    #
    #     layout = QVBoxLayout(self.Pymol_Widget)
    #     layout.addWidget(self.ProteinView)
    #     self.ProteinView.initial_pymol_visual()
    #     self.ProteinView.show()

    def start_pymol(self):
        # Creating the PyMolWidget
        try:
            self.ProteinView = PymolQtWidget(self)
            verticalLayoutProteinView = QVBoxLayout(self.Pymol_Widget)
            verticalLayoutProteinView.addWidget(self.ProteinView)
            self.setLayout(verticalLayoutProteinView)
            self.ProteinView.update()
            self.ProteinView.show()
            verticalLayoutProteinView.setContentsMargins(0, 0, 0, 0)
            # self.deleteLayout(verticalLayoutProteinView)
            self.ProteinView.initial_pymol_visual()
        except Exception as instance:
            if self.Pymol_Widget.isVisible():
                QMessageBox.critical(self, 'Visualize Problem!.', repr(instance) +
                                     "\n\nAn error occurred while loading the pdb file."
                                     "\n\nPlease select a pdb file for visualize by PyMol")

            else:
                QMessageBox.critical(self, 'PyMol can''be initialized', repr(instance) +
                                     "\n\nPyMol widget could not be initialized. "
                                     "Molecule viewer will not start, "
                                     "however the rest of the UI will be functional. "
                                     "Upgrading your graphics card drivers or reinstalling PyMol"
                                     "may solve this issue.")

    def deleteLayout(self, verticalLayoutProteinView):
        if verticalLayoutProteinView is not None:
            sip.delete(verticalLayoutProteinView)
    ## == > OPEN SOURCE PYMOL 2.4 INTEGRATION START
    ########################################################################
