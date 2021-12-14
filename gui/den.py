# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'MAIN_GUI.ui',
# licensing of 'MAIN_GUI.ui' applies.
#
# Created: Mon Dec 13 23:58:13 2021
#      by: pyside2-uic  running on PySide2 5.13.2
#
# WARNING! All changes made in this file will be lost!

from PySide2 import QtCore, QtGui, QtWidgets

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(1610, 1000)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.MinimumExpanding, QtWidgets.QSizePolicy.MinimumExpanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(MainWindow.sizePolicy().hasHeightForWidth())
        MainWindow.setSizePolicy(sizePolicy)
        MainWindow.setMinimumSize(QtCore.QSize(720, 1000))
        palette = QtGui.QPalette()
        brush = QtGui.QBrush(QtGui.QColor(255, 255, 255))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.WindowText, brush)
        brush = QtGui.QBrush(QtGui.QColor(0, 0, 0, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.Button, brush)
        brush = QtGui.QBrush(QtGui.QColor(66, 73, 90))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.Light, brush)
        brush = QtGui.QBrush(QtGui.QColor(55, 61, 75))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.Midlight, brush)
        brush = QtGui.QBrush(QtGui.QColor(22, 24, 30))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.Dark, brush)
        brush = QtGui.QBrush(QtGui.QColor(29, 32, 40))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.Mid, brush)
        brush = QtGui.QBrush(QtGui.QColor(210, 210, 210))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.Text, brush)
        brush = QtGui.QBrush(QtGui.QColor(255, 255, 255))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.BrightText, brush)
        brush = QtGui.QBrush(QtGui.QColor(255, 255, 255))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.ButtonText, brush)
        brush = QtGui.QBrush(QtGui.QColor(0, 0, 0, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.Base, brush)
        brush = QtGui.QBrush(QtGui.QColor(0, 0, 0, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.Window, brush)
        brush = QtGui.QBrush(QtGui.QColor(0, 0, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.Shadow, brush)
        brush = QtGui.QBrush(QtGui.QColor(85, 170, 255))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.Highlight, brush)
        brush = QtGui.QBrush(QtGui.QColor(85, 170, 255))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.Link, brush)
        brush = QtGui.QBrush(QtGui.QColor(255, 0, 127))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.LinkVisited, brush)
        brush = QtGui.QBrush(QtGui.QColor(22, 24, 30))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.AlternateBase, brush)
        brush = QtGui.QBrush(QtGui.QColor(44, 49, 60))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.ToolTipBase, brush)
        brush = QtGui.QBrush(QtGui.QColor(210, 210, 210))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.ToolTipText, brush)
        brush = QtGui.QBrush(QtGui.QColor(255, 255, 255))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Inactive, QtGui.QPalette.WindowText, brush)
        brush = QtGui.QBrush(QtGui.QColor(0, 0, 0, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Inactive, QtGui.QPalette.Button, brush)
        brush = QtGui.QBrush(QtGui.QColor(66, 73, 90))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Inactive, QtGui.QPalette.Light, brush)
        brush = QtGui.QBrush(QtGui.QColor(55, 61, 75))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Inactive, QtGui.QPalette.Midlight, brush)
        brush = QtGui.QBrush(QtGui.QColor(22, 24, 30))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Inactive, QtGui.QPalette.Dark, brush)
        brush = QtGui.QBrush(QtGui.QColor(29, 32, 40))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Inactive, QtGui.QPalette.Mid, brush)
        brush = QtGui.QBrush(QtGui.QColor(210, 210, 210))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Inactive, QtGui.QPalette.Text, brush)
        brush = QtGui.QBrush(QtGui.QColor(255, 255, 255))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Inactive, QtGui.QPalette.BrightText, brush)
        brush = QtGui.QBrush(QtGui.QColor(255, 255, 255))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Inactive, QtGui.QPalette.ButtonText, brush)
        brush = QtGui.QBrush(QtGui.QColor(0, 0, 0, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Inactive, QtGui.QPalette.Base, brush)
        brush = QtGui.QBrush(QtGui.QColor(0, 0, 0, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Inactive, QtGui.QPalette.Window, brush)
        brush = QtGui.QBrush(QtGui.QColor(0, 0, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Inactive, QtGui.QPalette.Shadow, brush)
        brush = QtGui.QBrush(QtGui.QColor(85, 170, 255))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Inactive, QtGui.QPalette.Highlight, brush)
        brush = QtGui.QBrush(QtGui.QColor(85, 170, 255))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Inactive, QtGui.QPalette.Link, brush)
        brush = QtGui.QBrush(QtGui.QColor(255, 0, 127))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Inactive, QtGui.QPalette.LinkVisited, brush)
        brush = QtGui.QBrush(QtGui.QColor(22, 24, 30))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Inactive, QtGui.QPalette.AlternateBase, brush)
        brush = QtGui.QBrush(QtGui.QColor(44, 49, 60))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Inactive, QtGui.QPalette.ToolTipBase, brush)
        brush = QtGui.QBrush(QtGui.QColor(210, 210, 210))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Inactive, QtGui.QPalette.ToolTipText, brush)
        brush = QtGui.QBrush(QtGui.QColor(22, 24, 30))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Disabled, QtGui.QPalette.WindowText, brush)
        brush = QtGui.QBrush(QtGui.QColor(0, 0, 0, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Disabled, QtGui.QPalette.Button, brush)
        brush = QtGui.QBrush(QtGui.QColor(66, 73, 90))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Disabled, QtGui.QPalette.Light, brush)
        brush = QtGui.QBrush(QtGui.QColor(55, 61, 75))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Disabled, QtGui.QPalette.Midlight, brush)
        brush = QtGui.QBrush(QtGui.QColor(22, 24, 30))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Disabled, QtGui.QPalette.Dark, brush)
        brush = QtGui.QBrush(QtGui.QColor(29, 32, 40))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Disabled, QtGui.QPalette.Mid, brush)
        brush = QtGui.QBrush(QtGui.QColor(22, 24, 30))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Disabled, QtGui.QPalette.Text, brush)
        brush = QtGui.QBrush(QtGui.QColor(255, 255, 255))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Disabled, QtGui.QPalette.BrightText, brush)
        brush = QtGui.QBrush(QtGui.QColor(22, 24, 30))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Disabled, QtGui.QPalette.ButtonText, brush)
        brush = QtGui.QBrush(QtGui.QColor(0, 0, 0, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Disabled, QtGui.QPalette.Base, brush)
        brush = QtGui.QBrush(QtGui.QColor(0, 0, 0, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Disabled, QtGui.QPalette.Window, brush)
        brush = QtGui.QBrush(QtGui.QColor(0, 0, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Disabled, QtGui.QPalette.Shadow, brush)
        brush = QtGui.QBrush(QtGui.QColor(51, 153, 255))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Disabled, QtGui.QPalette.Highlight, brush)
        brush = QtGui.QBrush(QtGui.QColor(85, 170, 255))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Disabled, QtGui.QPalette.Link, brush)
        brush = QtGui.QBrush(QtGui.QColor(255, 0, 127))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Disabled, QtGui.QPalette.LinkVisited, brush)
        brush = QtGui.QBrush(QtGui.QColor(44, 49, 60))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Disabled, QtGui.QPalette.AlternateBase, brush)
        brush = QtGui.QBrush(QtGui.QColor(44, 49, 60))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Disabled, QtGui.QPalette.ToolTipBase, brush)
        brush = QtGui.QBrush(QtGui.QColor(210, 210, 210))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Disabled, QtGui.QPalette.ToolTipText, brush)
        MainWindow.setPalette(palette)
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(10)
        MainWindow.setFont(font)
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap(":/big_icons/icons/big_icons/style_icon_48x48.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        MainWindow.setWindowIcon(icon)
        MainWindow.setStyleSheet("")
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setMinimumSize(QtCore.QSize(0, 0))
        self.centralwidget.setStyleSheet("\n"
"color: rgb(210, 210, 210);")
        self.centralwidget.setObjectName("centralwidget")
        self.horizontalLayout = QtWidgets.QHBoxLayout(self.centralwidget)
        self.horizontalLayout.setSpacing(0)
        self.horizontalLayout.setContentsMargins(0, 0, 0, 0)
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.frame_main = QtWidgets.QFrame(self.centralwidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.frame_main.sizePolicy().hasHeightForWidth())
        self.frame_main.setSizePolicy(sizePolicy)
        self.frame_main.setMinimumSize(QtCore.QSize(0, 0))
        self.frame_main.setStyleSheet("")
        self.frame_main.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_main.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_main.setObjectName("frame_main")
        self.verticalLayout = QtWidgets.QVBoxLayout(self.frame_main)
        self.verticalLayout.setSpacing(0)
        self.verticalLayout.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout.setObjectName("verticalLayout")
        self.frame_top = QtWidgets.QFrame(self.frame_main)
        self.frame_top.setMinimumSize(QtCore.QSize(0, 60))
        self.frame_top.setMaximumSize(QtCore.QSize(16777215, 60))
        self.frame_top.setStyleSheet("")
        self.frame_top.setFrameShape(QtWidgets.QFrame.NoFrame)
        self.frame_top.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_top.setObjectName("frame_top")
        self.horizontalLayout_3 = QtWidgets.QHBoxLayout(self.frame_top)
        self.horizontalLayout_3.setSpacing(0)
        self.horizontalLayout_3.setContentsMargins(0, 0, 0, 0)
        self.horizontalLayout_3.setObjectName("horizontalLayout_3")
        self.frame_toggle = QtWidgets.QFrame(self.frame_top)
        self.frame_toggle.setMaximumSize(QtCore.QSize(70, 16777215))
        self.frame_toggle.setStyleSheet("background-color: rgb(27, 29, 35);")
        self.frame_toggle.setFrameShape(QtWidgets.QFrame.NoFrame)
        self.frame_toggle.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_toggle.setObjectName("frame_toggle")
        self.verticalLayout_3 = QtWidgets.QVBoxLayout(self.frame_toggle)
        self.verticalLayout_3.setSpacing(0)
        self.verticalLayout_3.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout_3.setObjectName("verticalLayout_3")
        self.btn_toggle_menu = QtWidgets.QPushButton(self.frame_toggle)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.btn_toggle_menu.sizePolicy().hasHeightForWidth())
        self.btn_toggle_menu.setSizePolicy(sizePolicy)
        self.btn_toggle_menu.setStyleSheet("QPushButton {\n"
"    background-image: url(:/24x24/icons/24x24/cil-menu.png);\n"
"    background-position: center;\n"
"    background-repeat: no-reperat;\n"
"    border: none;\n"
"    background-color: rgb(27, 29, 35);\n"
"}\n"
"QPushButton:hover {\n"
"    background-color: rgb(33, 37, 43);\n"
"}\n"
"QPushButton:pressed {    \n"
"    background-color: rgb(85, 170, 255);\n"
"}")
        self.btn_toggle_menu.setText("")
        self.btn_toggle_menu.setObjectName("btn_toggle_menu")
        self.verticalLayout_3.addWidget(self.btn_toggle_menu)
        self.horizontalLayout_3.addWidget(self.frame_toggle)
        self.frame_top_right = QtWidgets.QFrame(self.frame_top)
        self.frame_top_right.setMinimumSize(QtCore.QSize(0, 45))
        self.frame_top_right.setMaximumSize(QtCore.QSize(16777215, 16777215))
        font = QtGui.QFont()
        font.setPointSize(8)
        self.frame_top_right.setFont(font)
        self.frame_top_right.setStyleSheet("background: transparent;")
        self.frame_top_right.setFrameShape(QtWidgets.QFrame.NoFrame)
        self.frame_top_right.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_top_right.setLineWidth(1)
        self.frame_top_right.setObjectName("frame_top_right")
        self.verticalLayout_2 = QtWidgets.QVBoxLayout(self.frame_top_right)
        self.verticalLayout_2.setSpacing(0)
        self.verticalLayout_2.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.frame_top_btns = QtWidgets.QFrame(self.frame_top_right)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.frame_top_btns.sizePolicy().hasHeightForWidth())
        self.frame_top_btns.setSizePolicy(sizePolicy)
        self.frame_top_btns.setMinimumSize(QtCore.QSize(0, 25))
        self.frame_top_btns.setMaximumSize(QtCore.QSize(16777215, 35))
        self.frame_top_btns.setStyleSheet("background-color: rgba(33, 37, 43, 150);\n"
"")
        self.frame_top_btns.setFrameShape(QtWidgets.QFrame.NoFrame)
        self.frame_top_btns.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_top_btns.setObjectName("frame_top_btns")
        self.horizontalLayout_4 = QtWidgets.QHBoxLayout(self.frame_top_btns)
        self.horizontalLayout_4.setSpacing(0)
        self.horizontalLayout_4.setContentsMargins(0, 0, 0, 0)
        self.horizontalLayout_4.setObjectName("horizontalLayout_4")
        self.frame_label_top_btns = QtWidgets.QFrame(self.frame_top_btns)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.frame_label_top_btns.sizePolicy().hasHeightForWidth())
        self.frame_label_top_btns.setSizePolicy(sizePolicy)
        self.frame_label_top_btns.setMinimumSize(QtCore.QSize(0, 0))
        self.frame_label_top_btns.setMaximumSize(QtCore.QSize(16777215, 16777215))
        self.frame_label_top_btns.setStyleSheet("border-color: rgb(0, 0, 0);\n"
"")
        self.frame_label_top_btns.setFrameShape(QtWidgets.QFrame.NoFrame)
        self.frame_label_top_btns.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_label_top_btns.setObjectName("frame_label_top_btns")
        self.horizontalLayout_10 = QtWidgets.QHBoxLayout(self.frame_label_top_btns)
        self.horizontalLayout_10.setSpacing(0)
        self.horizontalLayout_10.setContentsMargins(8, 0, 10, 0)
        self.horizontalLayout_10.setObjectName("horizontalLayout_10")
        self.frame_icon_top_bar = QtWidgets.QFrame(self.frame_label_top_btns)
        self.frame_icon_top_bar.setMaximumSize(QtCore.QSize(30, 23))
        self.frame_icon_top_bar.setStyleSheet("background: transparent;\n"
"background-image: url(:/16x16/icons/16x16/cil-terminal.png);\n"
"background-position: center;\n"
"background-repeat: no-repeat;")
        self.frame_icon_top_bar.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_icon_top_bar.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_icon_top_bar.setObjectName("frame_icon_top_bar")
        self.horizontalLayout_10.addWidget(self.frame_icon_top_bar)
        self.label_title_bar_top = QtWidgets.QLabel(self.frame_label_top_btns)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_title_bar_top.sizePolicy().hasHeightForWidth())
        self.label_title_bar_top.setSizePolicy(sizePolicy)
        self.label_title_bar_top.setMinimumSize(QtCore.QSize(0, 0))
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(10)
        font.setWeight(75)
        font.setBold(True)
        self.label_title_bar_top.setFont(font)
        self.label_title_bar_top.setStyleSheet("background-color: rgb(31, 36, 41);\n"
"\n"
"margin-left: 5px;\n"
"")
        self.label_title_bar_top.setObjectName("label_title_bar_top")
        self.horizontalLayout_10.addWidget(self.label_title_bar_top)
        self.horizontalLayout_4.addWidget(self.frame_label_top_btns)
        self.frame_btns_right = QtWidgets.QFrame(self.frame_top_btns)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.frame_btns_right.sizePolicy().hasHeightForWidth())
        self.frame_btns_right.setSizePolicy(sizePolicy)
        self.frame_btns_right.setMinimumSize(QtCore.QSize(0, 0))
        self.frame_btns_right.setMaximumSize(QtCore.QSize(120, 16777215))
        self.frame_btns_right.setFrameShape(QtWidgets.QFrame.NoFrame)
        self.frame_btns_right.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_btns_right.setObjectName("frame_btns_right")
        self.horizontalLayout_5 = QtWidgets.QHBoxLayout(self.frame_btns_right)
        self.horizontalLayout_5.setSpacing(0)
        self.horizontalLayout_5.setContentsMargins(0, 0, 0, 0)
        self.horizontalLayout_5.setObjectName("horizontalLayout_5")
        self.btn_minimize = QtWidgets.QPushButton(self.frame_btns_right)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.btn_minimize.sizePolicy().hasHeightForWidth())
        self.btn_minimize.setSizePolicy(sizePolicy)
        self.btn_minimize.setStyleSheet("QPushButton {    \n"
"    border: none;\n"
"    background-color: transparent;\n"
"}\n"
"QPushButton:hover {\n"
"    background-color: rgb(44, 49, 60)\n"
"}\n"
"QPushButton:pressed {    \n"
"    background-color: rgb(85, 170, 255);\n"
"}")
        self.btn_minimize.setText("")
        icon1 = QtGui.QIcon()
        icon1.addPixmap(QtGui.QPixmap(":/16x16/icons/16x16/cil-window-minimize.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.btn_minimize.setIcon(icon1)
        self.btn_minimize.setObjectName("btn_minimize")
        self.horizontalLayout_5.addWidget(self.btn_minimize)
        self.btn_maximize_restore = QtWidgets.QPushButton(self.frame_btns_right)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.btn_maximize_restore.sizePolicy().hasHeightForWidth())
        self.btn_maximize_restore.setSizePolicy(sizePolicy)
        self.btn_maximize_restore.setStyleSheet("QPushButton {    \n"
"    border: none;\n"
"    background-color: transparent;\n"
"}\n"
"QPushButton:hover {\n"
"    background-color: rgb(44, 49, 60)\n"
"}\n"
"QPushButton:pressed {    \n"
"    background-color: rgb(85, 170, 255);\n"
"}")
        self.btn_maximize_restore.setText("")
        icon2 = QtGui.QIcon()
        icon2.addPixmap(QtGui.QPixmap(":/16x16/icons/16x16/cil-window-maximize.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.btn_maximize_restore.setIcon(icon2)
        self.btn_maximize_restore.setObjectName("btn_maximize_restore")
        self.horizontalLayout_5.addWidget(self.btn_maximize_restore)
        self.btn_close = QtWidgets.QPushButton(self.frame_btns_right)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.btn_close.sizePolicy().hasHeightForWidth())
        self.btn_close.setSizePolicy(sizePolicy)
        self.btn_close.setStyleSheet("QPushButton {    \n"
"    border: none;\n"
"    background-color: transparent;\n"
"}\n"
"QPushButton:hover {\n"
"    background-color: rgb(44, 49, 60)\n"
"}\n"
"QPushButton:pressed {    \n"
"    background-color: rgb(85, 170, 255);\n"
"}")
        self.btn_close.setText("")
        icon3 = QtGui.QIcon()
        icon3.addPixmap(QtGui.QPixmap(":/16x16/icons/16x16/cil-x.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.btn_close.setIcon(icon3)
        self.btn_close.setObjectName("btn_close")
        self.horizontalLayout_5.addWidget(self.btn_close)
        self.horizontalLayout_4.addWidget(self.frame_btns_right)
        self.verticalLayout_2.addWidget(self.frame_top_btns)
        self.frame_top_info = QtWidgets.QFrame(self.frame_top_right)
        self.frame_top_info.setMinimumSize(QtCore.QSize(0, 25))
        self.frame_top_info.setMaximumSize(QtCore.QSize(16777215, 25))
        self.frame_top_info.setStyleSheet("/*background-color: rgb(39, 44, 54);*/\n"
"background-color: rgb(68, 72, 76);")
        self.frame_top_info.setFrameShape(QtWidgets.QFrame.NoFrame)
        self.frame_top_info.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_top_info.setObjectName("frame_top_info")
        self.horizontalLayout_8 = QtWidgets.QHBoxLayout(self.frame_top_info)
        self.horizontalLayout_8.setSpacing(0)
        self.horizontalLayout_8.setContentsMargins(10, 0, 10, 0)
        self.horizontalLayout_8.setObjectName("horizontalLayout_8")
        self.label_top_info_1 = QtWidgets.QLabel(self.frame_top_info)
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        self.label_top_info_1.setFont(font)
        self.label_top_info_1.setStyleSheet("color: rgb(98, 103, 111);")
        self.label_top_info_1.setObjectName("label_top_info_1")
        self.horizontalLayout_8.addWidget(self.label_top_info_1)
        self.label_top_info_2 = QtWidgets.QLabel(self.frame_top_info)
        self.label_top_info_2.setMinimumSize(QtCore.QSize(0, 0))
        self.label_top_info_2.setMaximumSize(QtCore.QSize(250, 16777215))
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setWeight(75)
        font.setBold(True)
        self.label_top_info_2.setFont(font)
        self.label_top_info_2.setStyleSheet("color: rgb(98, 103, 111);")
        self.label_top_info_2.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.label_top_info_2.setObjectName("label_top_info_2")
        self.horizontalLayout_8.addWidget(self.label_top_info_2)
        self.verticalLayout_2.addWidget(self.frame_top_info)
        self.horizontalLayout_3.addWidget(self.frame_top_right)
        self.verticalLayout.addWidget(self.frame_top)
        self.frame_center = QtWidgets.QFrame(self.frame_main)
        self.frame_center.setMinimumSize(QtCore.QSize(0, 0))
        self.frame_center.setStyleSheet("background-color: rgb(40, 44, 52);")
        self.frame_center.setFrameShape(QtWidgets.QFrame.NoFrame)
        self.frame_center.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_center.setObjectName("frame_center")
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout(self.frame_center)
        self.horizontalLayout_2.setSpacing(0)
        self.horizontalLayout_2.setContentsMargins(0, 0, 0, 0)
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.frame_left_menu = QtWidgets.QFrame(self.frame_center)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.frame_left_menu.sizePolicy().hasHeightForWidth())
        self.frame_left_menu.setSizePolicy(sizePolicy)
        self.frame_left_menu.setMinimumSize(QtCore.QSize(0, 0))
        self.frame_left_menu.setMaximumSize(QtCore.QSize(70, 16777215))
        self.frame_left_menu.setLayoutDirection(QtCore.Qt.LeftToRight)
        self.frame_left_menu.setStyleSheet("background-color: rgb(27, 29, 35);")
        self.frame_left_menu.setFrameShape(QtWidgets.QFrame.NoFrame)
        self.frame_left_menu.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_left_menu.setObjectName("frame_left_menu")
        self.verticalLayout_5 = QtWidgets.QVBoxLayout(self.frame_left_menu)
        self.verticalLayout_5.setSpacing(1)
        self.verticalLayout_5.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout_5.setObjectName("verticalLayout_5")
        self.frame_menus = QtWidgets.QFrame(self.frame_left_menu)
        self.frame_menus.setFrameShape(QtWidgets.QFrame.NoFrame)
        self.frame_menus.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_menus.setObjectName("frame_menus")
        self.layout_menus = QtWidgets.QVBoxLayout(self.frame_menus)
        self.layout_menus.setSpacing(0)
        self.layout_menus.setContentsMargins(0, 0, 0, 0)
        self.layout_menus.setObjectName("layout_menus")
        self.verticalLayout_5.addWidget(self.frame_menus)
        self.frame_extra_menus = QtWidgets.QFrame(self.frame_left_menu)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.frame_extra_menus.sizePolicy().hasHeightForWidth())
        self.frame_extra_menus.setSizePolicy(sizePolicy)
        self.frame_extra_menus.setFrameShape(QtWidgets.QFrame.NoFrame)
        self.frame_extra_menus.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_extra_menus.setObjectName("frame_extra_menus")
        self.layout_menu_bottom = QtWidgets.QVBoxLayout(self.frame_extra_menus)
        self.layout_menu_bottom.setSpacing(10)
        self.layout_menu_bottom.setContentsMargins(0, 0, 0, 25)
        self.layout_menu_bottom.setObjectName("layout_menu_bottom")
        self.label_user_icon = QtWidgets.QLabel(self.frame_extra_menus)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Maximum, QtWidgets.QSizePolicy.Maximum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_user_icon.sizePolicy().hasHeightForWidth())
        self.label_user_icon.setSizePolicy(sizePolicy)
        self.label_user_icon.setMinimumSize(QtCore.QSize(60, 60))
        self.label_user_icon.setMaximumSize(QtCore.QSize(60, 60))
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(12)
        self.label_user_icon.setFont(font)
        self.label_user_icon.setStyleSheet("QLabel {\n"
"    border-radius: 30px;\n"
"    background-color: rgb(44, 49, 60);\n"
"    border: 5px solid rgb(39, 44, 54);\n"
"}")
        self.label_user_icon.setAlignment(QtCore.Qt.AlignCenter)
        self.label_user_icon.setObjectName("label_user_icon")
        self.layout_menu_bottom.addWidget(self.label_user_icon)
        self.verticalLayout_5.addWidget(self.frame_extra_menus)
        self.horizontalLayout_2.addWidget(self.frame_left_menu)
        self.frame_content_right = QtWidgets.QFrame(self.frame_center)
        self.frame_content_right.setMinimumSize(QtCore.QSize(0, 0))
        self.frame_content_right.setStyleSheet("background-color: rgb(44, 49, 60);")
        self.frame_content_right.setFrameShape(QtWidgets.QFrame.NoFrame)
        self.frame_content_right.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_content_right.setObjectName("frame_content_right")
        self.verticalLayout_4 = QtWidgets.QVBoxLayout(self.frame_content_right)
        self.verticalLayout_4.setSpacing(0)
        self.verticalLayout_4.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout_4.setObjectName("verticalLayout_4")
        self.frame_content = QtWidgets.QFrame(self.frame_content_right)
        self.frame_content.setMinimumSize(QtCore.QSize(0, 820))
        self.frame_content.setMaximumSize(QtCore.QSize(16777215, 16777215))
        self.frame_content.setFrameShape(QtWidgets.QFrame.NoFrame)
        self.frame_content.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_content.setObjectName("frame_content")
        self.verticalLayout_9 = QtWidgets.QVBoxLayout(self.frame_content)
        self.verticalLayout_9.setSpacing(0)
        self.verticalLayout_9.setContentsMargins(5, 5, 5, 5)
        self.verticalLayout_9.setObjectName("verticalLayout_9")
        self.stackedWidget = QtWidgets.QStackedWidget(self.frame_content)
        self.stackedWidget.setMinimumSize(QtCore.QSize(0, 0))
        self.stackedWidget.setStyleSheet("background: transparent;")
        self.stackedWidget.setObjectName("stackedWidget")
        self.page_home = QtWidgets.QWidget()
        self.page_home.setObjectName("page_home")
        self.verticalLayout_10 = QtWidgets.QVBoxLayout(self.page_home)
        self.verticalLayout_10.setObjectName("verticalLayout_10")
        self.horizontalLayout_12 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_12.setObjectName("horizontalLayout_12")
        self.Perturbation_TabWidget = QtWidgets.QTabWidget(self.page_home)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.Perturbation_TabWidget.sizePolicy().hasHeightForWidth())
        self.Perturbation_TabWidget.setSizePolicy(sizePolicy)
        self.Perturbation_TabWidget.setMinimumSize(QtCore.QSize(600, 0))
        font = QtGui.QFont()
        font.setFamily("Segoe UI Semibold")
        font.setPointSize(8)
        self.Perturbation_TabWidget.setFont(font)
        self.Perturbation_TabWidget.setStyleSheet("\n"
"\n"
"QTabWidget::pane {\n"
"    border: 1px solid black;\n"
"    background: white;\n"
"    background-color: rgb(44, 49, 60);\n"
"    border-radius: 5px;\n"
"    border-top-color: rgb(157, 90, 198);\n"
"    border-right-color: rgb(157, 90, 198);\n"
"    border-bottom-color: rgb(157, 90, 198);\n"
"    border-left-color: rgb(157, 90, 198);\n"
"}\n"
"\n"
"QTabWidget::tab-bar:top {\n"
"    top: 0px;\n"
"\n"
"}\n"
"\n"
"QTabWidget::tab-bar:bottom {\n"
"    bottom: 0px;\n"
"}\n"
"\n"
"QTabWidget::tab-bar:left {\n"
"    right: 0px;\n"
"}\n"
"\n"
"QTabWidget::tab-bar:right {\n"
"    left: 0px;\n"
"}\n"
"\n"
"QTabBar::tab {\n"
"    border: 2px;\n"
"    \n"
"    border-color: rgb(0, 0, 0);\n"
"    border-radius: 2px;\n"
"    background-color: rgb(255, 17, 100);\n"
"\n"
"    \n"
"}\n"
"\n"
"\n"
"QTabBar::tab:selected {\n"
"    background: white;\n"
"    background-color: rgb(200, 5, 100);\n"
"    \n"
"    \n"
"}\n"
"\n"
"QTabBar::tab:!selected {\n"
"    background: white;\n"
"    \n"
"    background-color: rgb(141, 2, 72);\n"
"\n"
"    \n"
"}\n"
"\n"
"QTabBar::tab:!selected:hover {\n"
"    background: #999;\n"
"}\n"
"\n"
"QTabBar::tab:top:!selected {\n"
"    margin-top: 3px;\n"
"}\n"
"\n"
"QTabBar::tab:bottom:!selected {\n"
"    margin-bottom: 3px;\n"
"}\n"
"\n"
"QTabBar::tab:top, QTabBar::tab:bottom {\n"
"    min-width: 8ex;\n"
"    margin-right: -1px;\n"
"    padding: 5px 10px 5px 10px;\n"
"}\n"
"\n"
"QTabBar::tab:top:selected {\n"
"    border-bottom-color: none;\n"
"}\n"
"\n"
"QTabBar::tab:bottom:selected {\n"
"    border-top-color: none;\n"
"}\n"
"\n"
"QTabBar::tab:top:last, QTabBar::tab:bottom:last,\n"
"QTabBar::tab:top:only-one, QTabBar::tab:bottom:only-one {\n"
"    margin-right: 0;\n"
"}\n"
"\n"
"QTabBar::tab:left:!selected {\n"
"    margin-right: 3px;\n"
"}\n"
"\n"
"QTabBar::tab:right:!selected {\n"
"    margin-left: 3px;\n"
"}\n"
"\n"
"QTabBar::tab:left, QTabBar::tab:right {\n"
"    min-height: 8ex;\n"
"    margin-bottom: 1px;\n"
"    padding: 10px 5px 10px 5px;\n"
"}\n"
"\n"
"QTabBar::tab:left:selected {\n"
"    border-left-color: black;\n"
"}\n"
"\n"
"QTabBar::tab:right:selected {\n"
"    border-right-color: black;\n"
"}\n"
"\n"
"QTabBar::tab:left:last, QTabBar::tab:right:last,\n"
"QTabBar::tab:left:only-one, QTabBar::tab:right:only-one {\n"
"    margin-bottom: 0;\n"
"}")
        self.Perturbation_TabWidget.setObjectName("Perturbation_TabWidget")
        self.Perturb_Quick = QtWidgets.QWidget()
        self.Perturb_Quick.setObjectName("Perturb_Quick")
        self.verticalLayout_6 = QtWidgets.QVBoxLayout(self.Perturb_Quick)
        self.verticalLayout_6.setObjectName("verticalLayout_6")
        self.verticalLayout_7 = QtWidgets.QVBoxLayout()
        self.verticalLayout_7.setSpacing(10)
        self.verticalLayout_7.setSizeConstraint(QtWidgets.QLayout.SetDefaultConstraint)
        self.verticalLayout_7.setObjectName("verticalLayout_7")
        self.groupBox_14 = QtWidgets.QGroupBox(self.Perturb_Quick)
        self.groupBox_14.setMinimumSize(QtCore.QSize(0, 150))
        self.groupBox_14.setMaximumSize(QtCore.QSize(16777215, 141))
        font = QtGui.QFont()
        font.setFamily("Segoe UI Semibold")
        font.setPointSize(10)
        font.setUnderline(False)
        font.setStrikeOut(False)
        self.groupBox_14.setFont(font)
        self.groupBox_14.setStyleSheet("QGroupBox {\n"
"  \n"
"  border: 1px solid rgb(220, 220, 220);\n"
"  border-radius: 10px;\n"
"  padding: 4px;\n"
"  margin-top: 16px;\n"
"}\n"
"\n"
"QGroupBox::title {\n"
"  subcontrol-origin: margin;\n"
"  subcontrol-position: top left;\n"
"  left: 10px;\n"
"  padding-left: 3px;\n"
"  padding-right: 5px;\n"
"  padding-top: 8px;\n"
"  padding-bottom: 16px;\n"
"}\n"
"\n"
"QGroupBox::indicator {\n"
"  margin-left: 2px;\n"
"  height: 16px;\n"
"  width: 16px;\n"
"}")
        self.groupBox_14.setFlat(False)
        self.groupBox_14.setCheckable(False)
        self.groupBox_14.setObjectName("groupBox_14")
        self.horizontalLayout_42 = QtWidgets.QHBoxLayout(self.groupBox_14)
        self.horizontalLayout_42.setObjectName("horizontalLayout_42")
        self.gridLayout_2 = QtWidgets.QGridLayout()
        self.gridLayout_2.setSizeConstraint(QtWidgets.QLayout.SetFixedSize)
        self.gridLayout_2.setContentsMargins(3, 4, 0, 4)
        self.gridLayout_2.setObjectName("gridLayout_2")
        self.pdb_fix_clean = QtWidgets.QCheckBox(self.groupBox_14)
        font = QtGui.QFont()
        font.setPointSize(8)
        self.pdb_fix_clean.setFont(font)
        self.pdb_fix_clean.setStyleSheet("QCheckBox {\n"
"    spacing: 5px;\n"
"}\n"
"\n"
"QCheckBox::indicator {\n"
"    width: 15px;\n"
"    height: 15px;\n"
"}\n"
"\n"
"QCheckBox::indicator:unchecked {\n"
"border-image: url(:/16x16/icons/16x16/cil-circle.png);\n"
"\n"
"}\n"
"\n"
"\n"
"QCheckBox::indicator:checked {\n"
"    border-image: url(:/20x20/icons/20x20/cil-check.png);\n"
"}\n"
"\n"
"\n"
"")
        self.pdb_fix_clean.setChecked(True)
        self.pdb_fix_clean.setObjectName("pdb_fix_clean")
        self.gridLayout_2.addWidget(self.pdb_fix_clean, 1, 3, 1, 1)
        self.label_32 = QtWidgets.QLabel(self.groupBox_14)
        self.label_32.setMaximumSize(QtCore.QSize(83, 16777215))
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(10)
        self.label_32.setFont(font)
        self.label_32.setObjectName("label_32")
        self.gridLayout_2.addWidget(self.label_32, 2, 0, 1, 1)
        self.upload_pdb_Button = QtWidgets.QPushButton(self.groupBox_14)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.upload_pdb_Button.sizePolicy().hasHeightForWidth())
        self.upload_pdb_Button.setSizePolicy(sizePolicy)
        self.upload_pdb_Button.setMinimumSize(QtCore.QSize(106, 38))
        self.upload_pdb_Button.setMaximumSize(QtCore.QSize(106, 38))
        font = QtGui.QFont()
        font.setPointSize(10)
        font.setWeight(50)
        font.setBold(False)
        self.upload_pdb_Button.setFont(font)
        self.upload_pdb_Button.setStyleSheet("QPushButton {\n"
"    border: 2px solid rgb(52, 59, 72);\n"
"    border-radius: 5px;    \n"
"    background-color: rgb(0, 180, 95);\n"
"    margin-top: 2px;\n"
"    margin-bottom: 2px;\n"
"}\n"
"QPushButton:hover {\n"
"    background-color: rgb(0, 152, 84);\n"
"    border: 2px solid rgb(61, 70, 86);\n"
"}\n"
"QPushButton:pressed {    \n"
"    background-color: rgb(220, 5, 100);\n"
"    border: 2px solid rgb(43, 50, 61);\n"
"}\n"
"\n"
"QPushButton:disabled {\n"
"background-color: rgb(93, 83, 91);\n"
"}")
        self.upload_pdb_Button.setObjectName("upload_pdb_Button")
        self.gridLayout_2.addWidget(self.upload_pdb_Button, 2, 3, 1, 1)
        self.fetch_pdb_Button = QtWidgets.QPushButton(self.groupBox_14)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.fetch_pdb_Button.sizePolicy().hasHeightForWidth())
        self.fetch_pdb_Button.setSizePolicy(sizePolicy)
        self.fetch_pdb_Button.setMinimumSize(QtCore.QSize(0, 38))
        self.fetch_pdb_Button.setMaximumSize(QtCore.QSize(16777215, 38))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.fetch_pdb_Button.setFont(font)
        self.fetch_pdb_Button.setStyleSheet("QPushButton {\n"
"    border: 2px solid rgb(52, 59, 72);\n"
"    border-radius: 5px;    \n"
"    background-color: rgb(0, 180, 95);\n"
"    margin-top: 2px;\n"
"    margin-bottom: 2px;\n"
"}\n"
"QPushButton:hover {\n"
"    background-color: rgb(0, 152, 84);\n"
"    border: 2px solid rgb(61, 70, 86);\n"
"}\n"
"QPushButton:pressed {    \n"
"    background-color: rgb(220, 5, 100);\n"
"    border: 2px solid rgb(43, 50, 61);\n"
"}")
        self.fetch_pdb_Button.setObjectName("fetch_pdb_Button")
        self.gridLayout_2.addWidget(self.fetch_pdb_Button, 0, 3, 1, 1)
        self.label_2 = QtWidgets.QLabel(self.groupBox_14)
        font = QtGui.QFont()
        font.setPointSize(8)
        font.setWeight(50)
        font.setBold(False)
        self.label_2.setFont(font)
        self.label_2.setAlignment(QtCore.Qt.AlignCenter)
        self.label_2.setObjectName("label_2")
        self.gridLayout_2.addWidget(self.label_2, 1, 2, 1, 1)
        self.label = QtWidgets.QLabel(self.groupBox_14)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label.sizePolicy().hasHeightForWidth())
        self.label.setSizePolicy(sizePolicy)
        self.label.setMaximumSize(QtCore.QSize(83, 16777215))
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(10)
        font.setWeight(50)
        font.setBold(False)
        self.label.setFont(font)
        self.label.setObjectName("label")
        self.gridLayout_2.addWidget(self.label, 0, 0, 1, 1)
        self.upload_pdb_textEdit = QtWidgets.QTextEdit(self.groupBox_14)
        self.upload_pdb_textEdit.setMinimumSize(QtCore.QSize(375, 0))
        self.upload_pdb_textEdit.setMaximumSize(QtCore.QSize(382, 33))
        font = QtGui.QFont()
        font.setPointSize(8)
        self.upload_pdb_textEdit.setFont(font)
        self.upload_pdb_textEdit.setStyleSheet("QTextEdit {\n"
"    background-color: rgb(27, 29, 35);\n"
"    border-radius: 5px;\n"
"    padding: 10px;\n"
"}\n"
"QTextEdit:hover {\n"
"    border: 2px solid rgb(64, 71, 88);\n"
"}\n"
"QTextEdit:focus {\n"
"    border: 2px solid rgb(91, 101, 124);\n"
"}\n"
"\n"
"QTextEdit:disabled {\n"
"background-color: rgb(93, 83, 91);\n"
"}")
        self.upload_pdb_textEdit.setObjectName("upload_pdb_textEdit")
        self.gridLayout_2.addWidget(self.upload_pdb_textEdit, 2, 2, 1, 1)
        self.PDB_ID_lineEdit = QtWidgets.QLineEdit(self.groupBox_14)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.PDB_ID_lineEdit.sizePolicy().hasHeightForWidth())
        self.PDB_ID_lineEdit.setSizePolicy(sizePolicy)
        self.PDB_ID_lineEdit.setMinimumSize(QtCore.QSize(382, 33))
        self.PDB_ID_lineEdit.setMaximumSize(QtCore.QSize(382, 33))
        self.PDB_ID_lineEdit.setStyleSheet("QLineEdit {\n"
"    background-color: rgb(27, 29, 35);\n"
"    border-radius: 5px;\n"
"\n"
"}\n"
"QLineEdit:hover {\n"
"    border: 2px solid rgb(64, 71, 88);\n"
"}\n"
"QLineEdit:focus {\n"
"    border: 2px solid rgb(91, 101, 124);\n"
"}")
        self.PDB_ID_lineEdit.setObjectName("PDB_ID_lineEdit")
        self.gridLayout_2.addWidget(self.PDB_ID_lineEdit, 0, 2, 1, 1)
        self.horizontalLayout_42.addLayout(self.gridLayout_2)
        self.verticalLayout_7.addWidget(self.groupBox_14)
        self.groupBox_16 = QtWidgets.QGroupBox(self.Perturb_Quick)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.groupBox_16.sizePolicy().hasHeightForWidth())
        self.groupBox_16.setSizePolicy(sizePolicy)
        self.groupBox_16.setMinimumSize(QtCore.QSize(0, 90))
        self.groupBox_16.setMaximumSize(QtCore.QSize(16777215, 80))
        font = QtGui.QFont()
        font.setFamily("Segoe UI Semibold")
        font.setPointSize(10)
        self.groupBox_16.setFont(font)
        self.groupBox_16.setStyleSheet("QGroupBox {\n"
"  \n"
"  border: 1px solid rgb(220, 220, 220);\n"
"  border-radius: 10px;\n"
"  padding: 4px;\n"
"  margin-top: 16px;\n"
"}\n"
"\n"
"QGroupBox::title {\n"
"  subcontrol-origin: margin;\n"
"  subcontrol-position: top left;\n"
"  left: 10px;\n"
"  padding-left: 3px;\n"
"  padding-right: 5px;\n"
"  padding-top: 8px;\n"
"  padding-bottom: 16px;\n"
"}\n"
"\n"
"QGroupBox::indicator {\n"
"  margin-left: 2px;\n"
"  height: 16px;\n"
"  width: 16px;\n"
"}")
        self.groupBox_16.setObjectName("groupBox_16")
        self.verticalLayout_37 = QtWidgets.QVBoxLayout(self.groupBox_16)
        self.verticalLayout_37.setObjectName("verticalLayout_37")
        self.horizontalLayout_11 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_11.setContentsMargins(3, 4, 0, 4)
        self.horizontalLayout_11.setObjectName("horizontalLayout_11")
        self.label_31 = QtWidgets.QLabel(self.groupBox_16)
        self.label_31.setMinimumSize(QtCore.QSize(83, 0))
        self.label_31.setMaximumSize(QtCore.QSize(83, 16777215))
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(10)
        self.label_31.setFont(font)
        self.label_31.setObjectName("label_31")
        self.horizontalLayout_11.addWidget(self.label_31)
        self.Output_Folder_textEdit = QtWidgets.QTextEdit(self.groupBox_16)
        self.Output_Folder_textEdit.setMinimumSize(QtCore.QSize(382, 0))
        self.Output_Folder_textEdit.setMaximumSize(QtCore.QSize(382, 33))
        self.Output_Folder_textEdit.setStyleSheet("QTextEdit {\n"
"    background-color: rgb(27, 29, 35);\n"
"    border-radius: 5px;\n"
"    padding: 10px;\n"
"}\n"
"QTextEdit:hover {\n"
"    border: 2px solid rgb(64, 71, 88);\n"
"}\n"
"QTextEdit:focus {\n"
"    border: 2px solid rgb(91, 101, 124);\n"
"}")
        self.Output_Folder_textEdit.setObjectName("Output_Folder_textEdit")
        self.horizontalLayout_11.addWidget(self.Output_Folder_textEdit)
        self.Browse_Output_button = QtWidgets.QPushButton(self.groupBox_16)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.Browse_Output_button.sizePolicy().hasHeightForWidth())
        self.Browse_Output_button.setSizePolicy(sizePolicy)
        self.Browse_Output_button.setMinimumSize(QtCore.QSize(106, 38))
        self.Browse_Output_button.setMaximumSize(QtCore.QSize(106, 38))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.Browse_Output_button.setFont(font)
        self.Browse_Output_button.setStyleSheet("QPushButton {\n"
"    border: 2px solid rgb(52, 59, 72);\n"
"    border-radius: 5px;    \n"
"    background-color: rgb(0, 180, 95);\n"
"    margin-top:1px;\n"
"    margin-bottom: 1px;\n"
"\n"
"}\n"
"QPushButton:hover {\n"
"    background-color: rgb(0, 152, 84);\n"
"    border: 2px solid rgb(61, 70, 86);\n"
"}\n"
"QPushButton:pressed {    \n"
"    background-color: rgb(220, 5, 100);\n"
"    border: 2px solid rgb(43, 50, 61);\n"
"}")
        self.Browse_Output_button.setObjectName("Browse_Output_button")
        self.horizontalLayout_11.addWidget(self.Browse_Output_button)
        self.verticalLayout_37.addLayout(self.horizontalLayout_11)
        self.verticalLayout_7.addWidget(self.groupBox_16)
        self.groupBox_17 = QtWidgets.QGroupBox(self.Perturb_Quick)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.groupBox_17.sizePolicy().hasHeightForWidth())
        self.groupBox_17.setSizePolicy(sizePolicy)
        self.groupBox_17.setMinimumSize(QtCore.QSize(0, 120))
        self.groupBox_17.setMaximumSize(QtCore.QSize(16777215, 125))
        font = QtGui.QFont()
        font.setFamily("Segoe UI Semibold")
        font.setPointSize(10)
        self.groupBox_17.setFont(font)
        self.groupBox_17.setStyleSheet("QGroupBox {\n"
"  \n"
"  border: 1px solid rgb(220, 220, 220);\n"
"  border-radius: 10px;\n"
"  padding: 4px;\n"
"  margin-top: 16px;\n"
"}\n"
"\n"
"QGroupBox::title {\n"
"  subcontrol-origin: margin;\n"
"  subcontrol-position: top left;\n"
"  left: 10px;\n"
"  padding-left: 3px;\n"
"  padding-right: 5px;\n"
"  padding-top: 8px;\n"
"  padding-bottom: 16px;\n"
"}\n"
"\n"
"QGroupBox::indicator {\n"
"  margin-left: 2px;\n"
"  height: 16px;\n"
"  width: 16px;\n"
"}")
        self.groupBox_17.setObjectName("groupBox_17")
        self.horizontalLayout_40 = QtWidgets.QHBoxLayout(self.groupBox_17)
        self.horizontalLayout_40.setObjectName("horizontalLayout_40")
        self.gridLayout_4 = QtWidgets.QGridLayout()
        self.gridLayout_4.setContentsMargins(8, 3, 0, 3)
        self.gridLayout_4.setObjectName("gridLayout_4")
        self.label_35 = QtWidgets.QLabel(self.groupBox_17)
        self.label_35.setMinimumSize(QtCore.QSize(83, 0))
        self.label_35.setMaximumSize(QtCore.QSize(16777215, 16777215))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.label_35.setFont(font)
        self.label_35.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignVCenter)
        self.label_35.setObjectName("label_35")
        self.gridLayout_4.addWidget(self.label_35, 1, 0, 1, 1)
        self.All_CPU_checkBox = QtWidgets.QCheckBox(self.groupBox_17)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.All_CPU_checkBox.sizePolicy().hasHeightForWidth())
        self.All_CPU_checkBox.setSizePolicy(sizePolicy)
        self.All_CPU_checkBox.setMaximumSize(QtCore.QSize(200, 16777215))
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(10)
        self.All_CPU_checkBox.setFont(font)
        self.All_CPU_checkBox.setLayoutDirection(QtCore.Qt.LeftToRight)
        self.All_CPU_checkBox.setAutoFillBackground(False)
        self.All_CPU_checkBox.setStyleSheet("QCheckBox {\n"
"    spacing: 5px;\n"
"    margin-left: 15px;\n"
"}\n"
"\n"
"QCheckBox::indicator {\n"
"    width: 15px;\n"
"    height: 15px;\n"
"}\n"
"\n"
"QCheckBox::indicator:unchecked {\n"
"border-image: url(:/16x16/icons/16x16/cil-circle.png);\n"
"\n"
"}\n"
"\n"
"\n"
"QCheckBox::indicator:checked {\n"
"    border-image: url(:/20x20/icons/20x20/cil-check.png);\n"
"}\n"
"\n"
"\n"
"")
        self.All_CPU_checkBox.setObjectName("All_CPU_checkBox")
        self.gridLayout_4.addWidget(self.All_CPU_checkBox, 1, 2, 1, 1)
        self.long_simulation_time_unit = QtWidgets.QComboBox(self.groupBox_17)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.long_simulation_time_unit.sizePolicy().hasHeightForWidth())
        self.long_simulation_time_unit.setSizePolicy(sizePolicy)
        self.long_simulation_time_unit.setMinimumSize(QtCore.QSize(110, 0))
        self.long_simulation_time_unit.setMaximumSize(QtCore.QSize(16777215, 33))
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(10)
        self.long_simulation_time_unit.setFont(font)
        self.long_simulation_time_unit.setStyleSheet("QComboBox {\n"
"    background-color: rgb(27, 29, 35);\n"
"    border-radius: 5px;\n"
"    border: 2px solid rgb(27, 29, 35);\n"
"    padding: 1px 18px 1px 3px;\n"
"    padding-left: 10px;\n"
"    min-width: 6em;\n"
"}\n"
"\n"
"QComboBox:hover{\n"
"    border: 2px solid rgb(64, 71, 88);\n"
"}\n"
"\n"
"QComboBox QAbstractItemView{\n"
"    selection-background-color: rgb(200, 5, 100);\n"
"    selection-color: rgb(193, 125, 17);\n"
"}\n"
"\n"
"QComboBox:editable {\n"
"       selection-background-color: rgb(127, 5, 64);\n"
"    selection-color: rgb(255, 255, 255);\n"
"}\n"
"\n"
"QComboBox:on { /* shift the text when the popup opens */\n"
"    padding-top: 3px;\n"
"    padding-left: 4px;\n"
"    text-align: center;\n"
"}\n"
"\n"
"\n"
"\n"
"QComboBox:on { /* shift the text when the popup opens */\n"
"    padding-top: 2px;\n"
"    padding-left: 3px;\n"
"}\n"
"\n"
"QComboBox::drop-down {\n"
"    subcontrol-origin: padding;\n"
"    subcontrol-position: top right;\n"
"    width: 20px;\n"
"\n"
"    border-left-width: 2px;\n"
"    border-left-color: darkgray;\n"
"    border-left-style: solid; /* just a single line */\n"
"    border-top-right-radius: 3px; /* same radius as the QComboBox */\n"
"    border-bottom-right-radius: 3px;\n"
"}\n"
"\n"
"QComboBox::down-arrow {\n"
"image: url(:/16x16/icons/16x16/cil-arrow-bottom.png);\n"
"}\n"
"\n"
"QComboBox::down-arrow:on { /* shift the arrow when popup is open */\n"
"    top: 1px;\n"
"    left: 1px;\n"
"}\n"
"\n"
"")
        self.long_simulation_time_unit.setObjectName("long_simulation_time_unit")
        self.long_simulation_time_unit.addItem("")
        self.long_simulation_time_unit.addItem("")
        self.gridLayout_4.addWidget(self.long_simulation_time_unit, 0, 2, 1, 1)
        self.label_8 = QtWidgets.QLabel(self.groupBox_17)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_8.sizePolicy().hasHeightForWidth())
        self.label_8.setSizePolicy(sizePolicy)
        self.label_8.setMinimumSize(QtCore.QSize(83, 0))
        self.label_8.setMaximumSize(QtCore.QSize(16777215, 16777215))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.label_8.setFont(font)
        self.label_8.setStyleSheet("")
        self.label_8.setObjectName("label_8")
        self.gridLayout_4.addWidget(self.label_8, 0, 0, 1, 1)
        self.Number_CPU_spinBox_2 = QtWidgets.QSpinBox(self.groupBox_17)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.Number_CPU_spinBox_2.sizePolicy().hasHeightForWidth())
        self.Number_CPU_spinBox_2.setSizePolicy(sizePolicy)
        self.Number_CPU_spinBox_2.setMinimumSize(QtCore.QSize(171, 0))
        self.Number_CPU_spinBox_2.setMaximumSize(QtCore.QSize(16777215, 33))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.Number_CPU_spinBox_2.setFont(font)
        self.Number_CPU_spinBox_2.setLayoutDirection(QtCore.Qt.RightToLeft)
        self.Number_CPU_spinBox_2.setStyleSheet("QSpinBox{\n"
"    background-color: rgb(27, 29, 35);\n"
"    border-radius: 5px;\n"
"    border: 2px solid rgb(27, 29, 35);\n"
"    padding: 5px;\n"
"    padding-left: 5px;\n"
"}\n"
"QSpinBox:hover{\n"
"    border: 2px solid rgb(64, 71, 88);\n"
"}\n"
"QSpinBox QAbstractItemView {\n"
"    color: rgb(85, 170, 255);    \n"
"    background-color: rgb(27, 29, 35);\n"
"    padding: 5px;\n"
"    selection-background-color: rgb(39, 44, 54);\n"
"}")
        self.Number_CPU_spinBox_2.setMinimum(2)
        self.Number_CPU_spinBox_2.setObjectName("Number_CPU_spinBox_2")
        self.gridLayout_4.addWidget(self.Number_CPU_spinBox_2, 1, 1, 1, 1)
        self.run_duration_doubleSpinBox = QtWidgets.QDoubleSpinBox(self.groupBox_17)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.run_duration_doubleSpinBox.sizePolicy().hasHeightForWidth())
        self.run_duration_doubleSpinBox.setSizePolicy(sizePolicy)
        self.run_duration_doubleSpinBox.setMinimumSize(QtCore.QSize(0, 0))
        self.run_duration_doubleSpinBox.setMaximumSize(QtCore.QSize(16777215, 33))
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(10)
        self.run_duration_doubleSpinBox.setFont(font)
        self.run_duration_doubleSpinBox.setLayoutDirection(QtCore.Qt.RightToLeft)
        self.run_duration_doubleSpinBox.setStyleSheet("QDoubleSpinBox{\n"
"    background-color: rgb(27, 29, 35);\n"
"    border-radius: 5px;\n"
"    border: 2px solid rgb(27, 29, 35);\n"
"    padding: 5px;\n"
"    padding-left: 5px;\n"
"}\n"
"QDoubleSpinBox:hover{\n"
"    border: 2px solid rgb(64, 71, 88);\n"
"}\n"
"QDoubleSpinBox QAbstractItemView {\n"
"    color: rgb(85, 170, 255);    \n"
"    background-color: rgb(27, 29, 35);\n"
"    padding: 5px;\n"
"    selection-background-color: rgb(39, 44, 54);\n"
"}")
        self.run_duration_doubleSpinBox.setDecimals(3)
        self.run_duration_doubleSpinBox.setMinimum(0.0)
        self.run_duration_doubleSpinBox.setMaximum(500000.0)
        self.run_duration_doubleSpinBox.setProperty("value", 0.002)
        self.run_duration_doubleSpinBox.setObjectName("run_duration_doubleSpinBox")
        self.gridLayout_4.addWidget(self.run_duration_doubleSpinBox, 0, 1, 1, 1)
        self.horizontalLayout_40.addLayout(self.gridLayout_4)
        self.verticalLayout_7.addWidget(self.groupBox_17)
        self.groupBox_15 = QtWidgets.QGroupBox(self.Perturb_Quick)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.groupBox_15.sizePolicy().hasHeightForWidth())
        self.groupBox_15.setSizePolicy(sizePolicy)
        self.groupBox_15.setMinimumSize(QtCore.QSize(0, 210))
        self.groupBox_15.setMaximumSize(QtCore.QSize(16777215, 210))
        font = QtGui.QFont()
        font.setFamily("Segoe UI Semibold")
        font.setPointSize(10)
        self.groupBox_15.setFont(font)
        self.groupBox_15.setStyleSheet("QGroupBox {\n"
"  \n"
"  border: 1px solid rgb(220, 220, 220);\n"
"  border-radius: 10px;\n"
"  padding: 4px;\n"
"  margin-top: 16px;\n"
"}\n"
"\n"
"QGroupBox::title {\n"
"  subcontrol-origin: margin;\n"
"  subcontrol-position: top left;\n"
"  left: 10px;\n"
"  padding-left: 3px;\n"
"  padding-right: 5px;\n"
"  padding-top: 8px;\n"
"  padding-bottom: 16px;\n"
"}\n"
"\n"
"QGroupBox::indicator {\n"
"  margin-left: 2px;\n"
"  height: 16px;\n"
"  width: 16px;\n"
"}")
        self.groupBox_15.setObjectName("groupBox_15")
        self.horizontalLayout_39 = QtWidgets.QHBoxLayout(self.groupBox_15)
        self.horizontalLayout_39.setObjectName("horizontalLayout_39")
        self.gridLayout_3 = QtWidgets.QGridLayout()
        self.gridLayout_3.setContentsMargins(8, 4, 3, 4)
        self.gridLayout_3.setObjectName("gridLayout_3")
        self.label_5 = QtWidgets.QLabel(self.groupBox_15)
        self.label_5.setMinimumSize(QtCore.QSize(0, 0))
        self.label_5.setMaximumSize(QtCore.QSize(16777215, 16777215))
        font = QtGui.QFont()
        font.setPointSize(10)
        font.setWeight(50)
        font.setBold(False)
        self.label_5.setFont(font)
        self.label_5.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignVCenter)
        self.label_5.setObjectName("label_5")
        self.gridLayout_3.addWidget(self.label_5, 0, 0, 1, 1)
        self.res1_comboBox = QtWidgets.QComboBox(self.groupBox_15)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.res1_comboBox.sizePolicy().hasHeightForWidth())
        self.res1_comboBox.setSizePolicy(sizePolicy)
        self.res1_comboBox.setMinimumSize(QtCore.QSize(110, 33))
        self.res1_comboBox.setMaximumSize(QtCore.QSize(16777215, 33))
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(10)
        self.res1_comboBox.setFont(font)
        self.res1_comboBox.setStyleSheet("QComboBox {\n"
"    background-color: rgb(27, 29, 35);\n"
"    border-radius: 5px;\n"
"    border: 2px solid rgb(27, 29, 35);\n"
"    padding: 1px 18px 1px 3px;\n"
"    padding-left: 10px;\n"
"    min-width: 6em;\n"
"}\n"
"QComboBox:hover{\n"
"    border: 2px solid rgb(64, 71, 88);\n"
"}\n"
"\n"
"QComboBox QAbstractItemView{\n"
"    \n"
"    selection-background-color: rgb(200, 5, 100);\n"
"    \n"
"\n"
"}\n"
"\n"
"QComboBox:editable {\n"
"\n"
"    background: rgb(0, 180, 95);\n"
"}\n"
"\n"
"QComboBox:on { /* shift the text when the popup opens */\n"
"    padding-top: 3px;\n"
"    padding-left: 4px;\n"
"    text-align: center;\n"
"}\n"
"\n"
"\n"
"\n"
"QComboBox:on { /* shift the text when the popup opens */\n"
"    padding-top: 2px;\n"
"    padding-left: 3px;\n"
"}\n"
"\n"
"QComboBox::drop-down {\n"
"    subcontrol-origin: padding;\n"
"    subcontrol-position: top right;\n"
"    width: 20px;\n"
"\n"
"    border-left-width: 2px;\n"
"    border-left-color: darkgray;\n"
"    border-left-style: solid; /* just a single line */\n"
"    border-top-right-radius: 3px; /* same radius as the QComboBox */\n"
"    border-bottom-right-radius: 3px;\n"
"}\n"
"\n"
"QComboBox::down-arrow {\n"
"image: url(:/16x16/icons/16x16/cil-arrow-bottom.png);\n"
"}\n"
"\n"
"QComboBox::down-arrow:on { /* shift the arrow when popup is open */\n"
"    top: 1px;\n"
"    left: 1px;\n"
"}\n"
"\n"
"")
        self.res1_comboBox.setObjectName("res1_comboBox")
        self.gridLayout_3.addWidget(self.res1_comboBox, 0, 1, 1, 1)
        self.label_6 = QtWidgets.QLabel(self.groupBox_15)
        self.label_6.setMaximumSize(QtCore.QSize(16777215, 16777215))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.label_6.setFont(font)
        self.label_6.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignVCenter)
        self.label_6.setObjectName("label_6")
        self.gridLayout_3.addWidget(self.label_6, 3, 0, 1, 1)
        self.Number_CPU_spinBox = QtWidgets.QSpinBox(self.groupBox_15)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.Number_CPU_spinBox.sizePolicy().hasHeightForWidth())
        self.Number_CPU_spinBox.setSizePolicy(sizePolicy)
        self.Number_CPU_spinBox.setMinimumSize(QtCore.QSize(0, 33))
        self.Number_CPU_spinBox.setMaximumSize(QtCore.QSize(16777215, 33))
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(10)
        self.Number_CPU_spinBox.setFont(font)
        self.Number_CPU_spinBox.setLayoutDirection(QtCore.Qt.RightToLeft)
        self.Number_CPU_spinBox.setStyleSheet("QSpinBox{\n"
"    background-color: rgb(27, 29, 35);\n"
"    border-radius: 5px;\n"
"    border: 2px solid rgb(27, 29, 35);\n"
"    padding: 5px;\n"
"    padding-left: 5px;\n"
"}\n"
"QSpinBox:hover{\n"
"    border: 2px solid rgb(64, 71, 88);\n"
"}\n"
"QSpinBox QAbstractItemView {\n"
"    color: rgb(85, 170, 255);    \n"
"    background-color: rgb(27, 29, 35);\n"
"    padding: 5px;\n"
"    selection-background-color: rgb(39, 44, 54);\n"
"}")
        self.Number_CPU_spinBox.setMinimum(2)
        self.Number_CPU_spinBox.setObjectName("Number_CPU_spinBox")
        self.gridLayout_3.addWidget(self.Number_CPU_spinBox, 3, 1, 1, 1)
        self.checkBox_2 = QtWidgets.QCheckBox(self.groupBox_15)
        font = QtGui.QFont()
        font.setPointSize(10)
        self.checkBox_2.setFont(font)
        self.checkBox_2.setStyleSheet("QCheckBox {\n"
"    spacing: 5px;\n"
"    margin-left: 15px;\n"
"}\n"
"\n"
"QCheckBox::indicator {\n"
"    width: 15px;\n"
"    height: 15px;\n"
"}\n"
"\n"
"QCheckBox::indicator:unchecked {\n"
"border-image: url(:/16x16/icons/16x16/cil-circle.png);\n"
"\n"
"}\n"
"\n"
"\n"
"QCheckBox::indicator:checked {\n"
"    border-image: url(:/20x20/icons/20x20/cil-check.png);\n"
"}")
        self.checkBox_2.setObjectName("checkBox_2")
        self.gridLayout_3.addWidget(self.checkBox_2, 3, 2, 1, 1)
        self.horizontalLayout_13 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_13.setSpacing(15)
        self.horizontalLayout_13.setContentsMargins(2, -1, 1, -1)
        self.horizontalLayout_13.setObjectName("horizontalLayout_13")
        self.discard_residue_pushButton = QtWidgets.QPushButton(self.groupBox_15)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.discard_residue_pushButton.sizePolicy().hasHeightForWidth())
        self.discard_residue_pushButton.setSizePolicy(sizePolicy)
        self.discard_residue_pushButton.setMinimumSize(QtCore.QSize(50, 33))
        self.discard_residue_pushButton.setMaximumSize(QtCore.QSize(60, 33))
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(9)
        self.discard_residue_pushButton.setFont(font)
        self.discard_residue_pushButton.setStyleSheet("QPushButton {\n"
"    border: 2px solid rgb(52, 59, 72);\n"
"    border-radius: 5px;    \n"
"    background-color: rgb(0, 180, 95);\n"
"}\n"
"QPushButton:hover {\n"
"    background-color: rgb(0, 152, 84);\n"
"    border: 2px solid rgb(61, 70, 86);\n"
"}\n"
"QPushButton:pressed {    \n"
"    background-color: rgb(220, 5, 100);\n"
"    border: 2px solid rgb(43, 50, 61);\n"
"}")
        self.discard_residue_pushButton.setText("")
        icon4 = QtGui.QIcon()
        icon4.addPixmap(QtGui.QPixmap(":/16x16/icons/16x16/cil-chevron-double-left.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.discard_residue_pushButton.setIcon(icon4)
        self.discard_residue_pushButton.setObjectName("discard_residue_pushButton")
        self.horizontalLayout_13.addWidget(self.discard_residue_pushButton)
        self.add_residue_pushButton = QtWidgets.QPushButton(self.groupBox_15)
        self.add_residue_pushButton.setMinimumSize(QtCore.QSize(50, 33))
        self.add_residue_pushButton.setMaximumSize(QtCore.QSize(60, 33))
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(9)
        self.add_residue_pushButton.setFont(font)
        self.add_residue_pushButton.setStyleSheet("QPushButton {\n"
"    border: 2px solid rgb(52, 59, 72);\n"
"    border-radius: 5px;    \n"
"    background-color: rgb(0, 180, 95);\n"
"}\n"
"QPushButton:hover {\n"
"    background-color: rgb(0, 152, 84);\n"
"    border: 2px solid rgb(61, 70, 86);\n"
"}\n"
"QPushButton:pressed {    \n"
"    background-color: rgb(220, 5, 100);\n"
"    border: 2px solid rgb(43, 50, 61);\n"
"}")
        self.add_residue_pushButton.setText("")
        icon5 = QtGui.QIcon()
        icon5.addPixmap(QtGui.QPixmap(":/16x16/icons/16x16/cil-chevron-double-right.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.add_residue_pushButton.setIcon(icon5)
        self.add_residue_pushButton.setObjectName("add_residue_pushButton")
        self.horizontalLayout_13.addWidget(self.add_residue_pushButton)
        self.gridLayout_3.addLayout(self.horizontalLayout_13, 0, 2, 1, 1)
        self.run_duration_spinBox = QtWidgets.QSpinBox(self.groupBox_15)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.run_duration_spinBox.sizePolicy().hasHeightForWidth())
        self.run_duration_spinBox.setSizePolicy(sizePolicy)
        self.run_duration_spinBox.setMinimumSize(QtCore.QSize(170, 33))
        self.run_duration_spinBox.setMaximumSize(QtCore.QSize(16777215, 33))
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(10)
        self.run_duration_spinBox.setFont(font)
        self.run_duration_spinBox.setLayoutDirection(QtCore.Qt.RightToLeft)
        self.run_duration_spinBox.setStyleSheet("QSpinBox{\n"
"    background-color: rgb(27, 29, 35);\n"
"    border-radius: 5px;\n"
"    border: 2px solid rgb(27, 29, 35);\n"
"    padding: 5px;\n"
"    padding-left: 5px;\n"
"}\n"
"QSpinBox:hover{\n"
"    border: 2px solid rgb(64, 71, 88);\n"
"}\n"
"\n"
"QSpinBox QAbstractItemView {\n"
"    color: rgb(85, 170, 255);    \n"
"    background-color: rgb(27, 29, 35);\n"
"    padding: 5px;\n"
"    selection-background-color: rgb(39, 44, 54);\n"
"}")
        self.run_duration_spinBox.setMaximum(10000000)
        self.run_duration_spinBox.setProperty("value", 5)
        self.run_duration_spinBox.setObjectName("run_duration_spinBox")
        self.gridLayout_3.addWidget(self.run_duration_spinBox, 2, 1, 1, 1)
        self.perturb_time_unit_comboBox = QtWidgets.QComboBox(self.groupBox_15)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.perturb_time_unit_comboBox.sizePolicy().hasHeightForWidth())
        self.perturb_time_unit_comboBox.setSizePolicy(sizePolicy)
        self.perturb_time_unit_comboBox.setMinimumSize(QtCore.QSize(110, 33))
        self.perturb_time_unit_comboBox.setMaximumSize(QtCore.QSize(16777215, 33))
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(10)
        self.perturb_time_unit_comboBox.setFont(font)
        self.perturb_time_unit_comboBox.setStyleSheet("QComboBox {\n"
"    background-color: rgb(27, 29, 35);\n"
"    border-radius: 5px;\n"
"    border: 2px solid rgb(27, 29, 35);\n"
"    padding: 1px 18px 1px 3px;\n"
"    padding-left: 10px;\n"
"    min-width: 6em;\n"
"}\n"
"QComboBox:hover{\n"
"    border: 2px solid rgb(64, 71, 88);\n"
"}\n"
"\n"
"QComboBox QAbstractItemView{\n"
"    \n"
"    selection-background-color: rgb(200, 5, 100);\n"
"    \n"
"\n"
"}\n"
"\n"
"QComboBox:editable {\n"
"\n"
"    background: rgb(0, 180, 95);\n"
"}\n"
"\n"
"QComboBox:on { /* shift the text when the popup opens */\n"
"    padding-top: 3px;\n"
"    padding-left: 4px;\n"
"    text-align: center;\n"
"}\n"
"\n"
"\n"
"\n"
"QComboBox:on { /* shift the text when the popup opens */\n"
"    padding-top: 2px;\n"
"    padding-left: 3px;\n"
"}\n"
"\n"
"QComboBox::drop-down {\n"
"    subcontrol-origin: padding;\n"
"    subcontrol-position: top right;\n"
"    width: 20px;\n"
"\n"
"    border-left-width: 2px;\n"
"    border-left-color: darkgray;\n"
"    border-left-style: solid; /* just a single line */\n"
"    border-top-right-radius: 3px; /* same radius as the QComboBox */\n"
"    border-bottom-right-radius: 3px;\n"
"}\n"
"\n"
"QComboBox::down-arrow {\n"
"image: url(:/16x16/icons/16x16/cil-arrow-bottom.png);\n"
"}\n"
"\n"
"QComboBox::down-arrow:on { /* shift the arrow when popup is open */\n"
"    top: 1px;\n"
"    left: 1px;\n"
"}\n"
"\n"
"")
        self.perturb_time_unit_comboBox.setObjectName("perturb_time_unit_comboBox")
        self.perturb_time_unit_comboBox.addItem("")
        self.gridLayout_3.addWidget(self.perturb_time_unit_comboBox, 2, 2, 1, 1)
        self.label_4 = QtWidgets.QLabel(self.groupBox_15)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_4.sizePolicy().hasHeightForWidth())
        self.label_4.setSizePolicy(sizePolicy)
        self.label_4.setMaximumSize(QtCore.QSize(16777215, 16777215))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.label_4.setFont(font)
        self.label_4.setStyleSheet("QLabel{\n"
"margin-right: 20px;\n"
"}")
        self.label_4.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignVCenter)
        self.label_4.setObjectName("label_4")
        self.gridLayout_3.addWidget(self.label_4, 1, 0, 1, 1)
        self.label_3 = QtWidgets.QLabel(self.groupBox_15)
        self.label_3.setMaximumSize(QtCore.QSize(16777215, 16777215))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.label_3.setFont(font)
        self.label_3.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignVCenter)
        self.label_3.setObjectName("label_3")
        self.gridLayout_3.addWidget(self.label_3, 2, 0, 1, 1)
        self.R_factor_lineEdit = QtWidgets.QLineEdit(self.groupBox_15)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.R_factor_lineEdit.sizePolicy().hasHeightForWidth())
        self.R_factor_lineEdit.setSizePolicy(sizePolicy)
        self.R_factor_lineEdit.setMinimumSize(QtCore.QSize(170, 33))
        self.R_factor_lineEdit.setMaximumSize(QtCore.QSize(16777215, 33))
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(10)
        self.R_factor_lineEdit.setFont(font)
        self.R_factor_lineEdit.setStyleSheet("QLineEdit {\n"
"    background-color: rgb(27, 29, 35);\n"
"    border-radius: 5px;\n"
"    padding-left: 5px;\n"
"\n"
"}\n"
"QLineEdit:hover {\n"
"    border: 2px solid rgb(64, 71, 88);\n"
"}\n"
"QLineEdit:focus {\n"
"    border: 2px solid rgb(91, 101, 124);\n"
"}")
        self.R_factor_lineEdit.setObjectName("R_factor_lineEdit")
        self.gridLayout_3.addWidget(self.R_factor_lineEdit, 1, 1, 1, 1)
        self.horizontalLayout_39.addLayout(self.gridLayout_3)
        self.verticalLayout_17 = QtWidgets.QVBoxLayout()
        self.verticalLayout_17.setSpacing(8)
        self.verticalLayout_17.setContentsMargins(3, 4, 3, 8)
        self.verticalLayout_17.setObjectName("verticalLayout_17")
        self.label_7 = QtWidgets.QLabel(self.groupBox_15)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_7.sizePolicy().hasHeightForWidth())
        self.label_7.setSizePolicy(sizePolicy)
        self.label_7.setMaximumSize(QtCore.QSize(180, 16777215))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.label_7.setFont(font)
        self.label_7.setStyleSheet("QLabel {\n"
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
        self.label_7.setObjectName("label_7")
        self.verticalLayout_17.addWidget(self.label_7)
        self.selected_residues_listWidget = QtWidgets.QListWidget(self.groupBox_15)
        self.selected_residues_listWidget.setMinimumSize(QtCore.QSize(0, 125))
        self.selected_residues_listWidget.setMaximumSize(QtCore.QSize(180, 16777215))
        self.selected_residues_listWidget.setStyleSheet("QListWidget {\n"
"    background-color: rgb(27, 29, 35);\n"
"    border-radius: 5px;\n"
"    border: 2px solid rgb(27, 29, 35);\n"
"    padding: 1px 18px 1px 3px;\n"
"    padding-left: 10px;\n"
"    selection-background-color: rgb(127, 5, 64);\n"
"    selection-color: rgb(255, 255, 255);\n"
"    font-size: 13px;\n"
"}\n"
"\n"
"QListWidget:hover{\n"
"    border: 2px solid rgb(64, 71, 88);\n"
"    selection-color: rgb(127, 5, 64);\n"
"\n"
"}\n"
"\n"
" QScrollBar:vertical\n"
"        {\n"
"            background-color: #323232;\n"
"            width: 15px;\n"
"            margin: 15px 3px 15px 3px;\n"
"            border: 1px transparent #2A2929;\n"
"            margin-right: 5px;\n"
"        }\n"
"        \n"
"        QScrollBar::handle:vertical\n"
"        {\n"
"            background-color: rgb(255, 17, 100);\n"
"            min-height: 5px;\n"
"            border-radius: 3px;\n"
"            \n"
"        }\n"
"        \n"
"        QScrollBar::sub-line:vertical\n"
"        {\n"
"            border: none;\n"
"            background: none;\n"
"            color: none;\n"
"        }\n"
"        \n"
"        QScrollBar::add-line:vertical\n"
"        {\n"
"            border: none;\n"
"            background: none;\n"
"            color: none;\n"
"        }\n"
"        \n"
"        QScrollBar::add-page:vertical, QScrollBar::sub-page:vertical\n"
"        {\n"
"            background: none;\n"
"        }\n"
"\n"
"")
        self.selected_residues_listWidget.setObjectName("selected_residues_listWidget")
        self.verticalLayout_17.addWidget(self.selected_residues_listWidget)
        self.horizontalLayout_39.addLayout(self.verticalLayout_17)
        self.verticalLayout_7.addWidget(self.groupBox_15)
        spacerItem = QtWidgets.QSpacerItem(20, 60, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum)
        self.verticalLayout_7.addItem(spacerItem)
        self.groupBox = QtWidgets.QGroupBox(self.Perturb_Quick)
        self.groupBox.setMinimumSize(QtCore.QSize(0, 60))
        self.groupBox.setMaximumSize(QtCore.QSize(16777215, 60))
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        self.groupBox.setFont(font)
        self.groupBox.setStyleSheet("")
        self.groupBox.setTitle("")
        self.groupBox.setObjectName("groupBox")
        self.horizontalLayout_43 = QtWidgets.QHBoxLayout(self.groupBox)
        self.horizontalLayout_43.setObjectName("horizontalLayout_43")
        self.horizontalLayout_16 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_16.setObjectName("horizontalLayout_16")
        self.quit_pushButton = QtWidgets.QPushButton(self.groupBox)
        self.quit_pushButton.setMinimumSize(QtCore.QSize(0, 33))
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(9)
        self.quit_pushButton.setFont(font)
        self.quit_pushButton.setStyleSheet("QPushButton {\n"
"    border: 2px solid rgb(52, 59, 72);\n"
"    border-radius: 5px;    \n"
"    background-color: rgb(0, 180, 95);\n"
"}\n"
"QPushButton:hover {\n"
"    background-color: rgb(0, 152, 84);\n"
"    border: 2px solid rgb(61, 70, 86);\n"
"}\n"
"QPushButton:pressed {    \n"
"    background-color: rgb(220, 5, 100);\n"
"    border: 2px solid rgb(43, 50, 61);\n"
"}")
        icon6 = QtGui.QIcon()
        icon6.addPixmap(QtGui.QPixmap(":/16x16/icons/16x16/cil-account-logout.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.quit_pushButton.setIcon(icon6)
        self.quit_pushButton.setObjectName("quit_pushButton")
        self.horizontalLayout_16.addWidget(self.quit_pushButton)
        self.stop_pushButton = QtWidgets.QPushButton(self.groupBox)
        self.stop_pushButton.setMinimumSize(QtCore.QSize(0, 33))
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(9)
        self.stop_pushButton.setFont(font)
        self.stop_pushButton.setStyleSheet("QPushButton {\n"
"    border: 2px solid rgb(52, 59, 72);\n"
"    border-radius: 5px;    \n"
"    background-color: rgb(0, 180, 95);\n"
"}\n"
"QPushButton:hover {\n"
"    background-color: rgb(0, 152, 84);\n"
"    border: 2px solid rgb(61, 70, 86);\n"
"}\n"
"QPushButton:pressed {    \n"
"    background-color: rgb(220, 5, 100);\n"
"    border: 2px solid rgb(43, 50, 61);\n"
"}")
        icon7 = QtGui.QIcon()
        icon7.addPixmap(QtGui.QPixmap(":/16x16/icons/16x16/cil-media-stop.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.stop_pushButton.setIcon(icon7)
        self.stop_pushButton.setObjectName("stop_pushButton")
        self.horizontalLayout_16.addWidget(self.stop_pushButton)
        self.pushButton_2 = QtWidgets.QPushButton(self.groupBox)
        self.pushButton_2.setMinimumSize(QtCore.QSize(0, 33))
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(9)
        self.pushButton_2.setFont(font)
        self.pushButton_2.setStyleSheet("QPushButton {\n"
"    border: 2px solid rgb(52, 59, 72);\n"
"    border-radius: 5px;    \n"
"    background-color: rgb(0, 180, 95);\n"
"}\n"
"QPushButton:hover {\n"
"    background-color: rgb(0, 152, 84);\n"
"    border: 2px solid rgb(61, 70, 86);\n"
"}\n"
"QPushButton:pressed {    \n"
"    background-color: rgb(220, 5, 100);\n"
"    border: 2px solid rgb(43, 50, 61);\n"
"}")
        icon8 = QtGui.QIcon()
        icon8.addPixmap(QtGui.QPixmap(":/16x16/icons/16x16/cil-save.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.pushButton_2.setIcon(icon8)
        self.pushButton_2.setObjectName("pushButton_2")
        self.horizontalLayout_16.addWidget(self.pushButton_2)
        self.Run = QtWidgets.QPushButton(self.groupBox)
        self.Run.setMinimumSize(QtCore.QSize(0, 33))
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(9)
        self.Run.setFont(font)
        self.Run.setStyleSheet("QPushButton {\n"
"    border: 2px solid rgb(52, 59, 72);\n"
"    border-radius: 5px;    \n"
"    background-color: rgb(0, 180, 95);\n"
"}\n"
"QPushButton:hover {\n"
"    background-color: rgb(255, 170, 0);\n"
"    border: 2px solid rgb(61, 70, 86);\n"
"}\n"
"QPushButton:pressed {    \n"
"    background-color: rgb(220, 5, 100);\n"
"    border: 2px solid rgb(43, 50, 61);\n"
"}\n"
"\n"
"")
        icon9 = QtGui.QIcon()
        icon9.addPixmap(QtGui.QPixmap(":/16x16/icons/16x16/cil-media-play.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.Run.setIcon(icon9)
        self.Run.setObjectName("Run")
        self.horizontalLayout_16.addWidget(self.Run)
        self.horizontalLayout_43.addLayout(self.horizontalLayout_16)
        self.verticalLayout_7.addWidget(self.groupBox)
        self.verticalLayout_6.addLayout(self.verticalLayout_7)
        self.Perturbation_TabWidget.addTab(self.Perturb_Quick, "")
        self.Perturb_Advanced = QtWidgets.QWidget()
        self.Perturb_Advanced.setObjectName("Perturb_Advanced")
        self.verticalLayout_21 = QtWidgets.QVBoxLayout(self.Perturb_Advanced)
        self.verticalLayout_21.setObjectName("verticalLayout_21")
        self.tabWidget = QtWidgets.QTabWidget(self.Perturb_Advanced)
        font = QtGui.QFont()
        font.setFamily("Segoe UI Semibold")
        font.setPointSize(8)
        self.tabWidget.setFont(font)
        self.tabWidget.setStyleSheet("QTabWidget::pane {\n"
"    border: 1px solid black;\n"
"    background: white;\n"
"    background-color: rgb(44, 49, 60);\n"
"    border-radius: 10px;\n"
"    border-top-color: rgb(157, 90, 198);\n"
"    border-right-color: rgb(0, 0, 0);\n"
"    border-bottom-color: none;\n"
"    border-left-color: rgb(0, 0, 0);\n"
"}\n"
"\n"
"QTabWidget::tab-bar:top {\n"
"    top: 0px;\n"
"\n"
"}\n"
"\n"
"QTabWidget::tab-bar:bottom {\n"
"    bottom: 0px;\n"
"}\n"
"\n"
"QTabWidget::tab-bar:left {\n"
"    right: 0px;\n"
"}\n"
"\n"
"QTabWidget::tab-bar:right {\n"
"    left: 0px;\n"
"}\n"
"\n"
"QTabBar::tab {\n"
"    border: 2px;\n"
"    \n"
"    border-color: rgb(0, 0, 0);\n"
"    border-radius: 2px;\n"
"    background-color: rgb(255, 17, 100);\n"
"\n"
"    \n"
"}\n"
"\n"
"\n"
"QTabBar::tab:selected {\n"
"    background: white;\n"
"    background-color: rgb(200, 5, 100);\n"
"    \n"
"    \n"
"}\n"
"\n"
"QTabBar::tab:!selected {\n"
"    background: white;\n"
"    \n"
"    background-color: rgb(141, 2, 72);\n"
"\n"
"    \n"
"}\n"
"\n"
"QTabBar::tab:!selected:hover {\n"
"    background: #999;\n"
"}\n"
"\n"
"QTabBar::tab:top:!selected {\n"
"    margin-top: 3px;\n"
"}\n"
"\n"
"QTabBar::tab:bottom:!selected {\n"
"    margin-bottom: 3px;\n"
"}\n"
"\n"
"QTabBar::tab:top, QTabBar::tab:bottom {\n"
"    min-width: 8ex;\n"
"    margin-right: -1px;\n"
"    padding: 5px 10px 5px 10px;\n"
"}\n"
"\n"
"QTabBar::tab:top:selected {\n"
"    border-bottom-color: none;\n"
"}\n"
"\n"
"QTabBar::tab:bottom:selected {\n"
"    border-top-color: none;\n"
"}\n"
"\n"
"QTabBar::tab:top:last, QTabBar::tab:bottom:last,\n"
"QTabBar::tab:top:only-one, QTabBar::tab:bottom:only-one {\n"
"    margin-right: 0;\n"
"}\n"
"\n"
"QTabBar::tab:left:!selected {\n"
"    margin-right: 3px;\n"
"}\n"
"\n"
"QTabBar::tab:right:!selected {\n"
"    margin-left: 3px;\n"
"}\n"
"\n"
"QTabBar::tab:left, QTabBar::tab:right {\n"
"    min-height: 8ex;\n"
"    margin-bottom: 1px;\n"
"    padding: 5px 5px 5px 5px;\n"
"}\n"
"\n"
"QTabBar::tab:left:selected {\n"
"    border-left-color: black;\n"
"}\n"
"\n"
"QTabBar::tab:right:selected {\n"
"    border-right-color: black;\n"
"}\n"
"\n"
"QTabBar::tab:left:last, QTabBar::tab:right:last,\n"
"QTabBar::tab:left:only-one, QTabBar::tab:right:only-one {\n"
"    margin-bottom: 0;\n"
"}")
        self.tabWidget.setTabPosition(QtWidgets.QTabWidget.North)
        self.tabWidget.setTabShape(QtWidgets.QTabWidget.Rounded)
        self.tabWidget.setMovable(True)
        self.tabWidget.setObjectName("tabWidget")
        self.tab = QtWidgets.QWidget()
        self.tab.setObjectName("tab")
        self.layoutWidget = QtWidgets.QWidget(self.tab)
        self.layoutWidget.setGeometry(QtCore.QRect(10, 10, 581, 261))
        self.layoutWidget.setObjectName("layoutWidget")
        self.verticalLayout_8 = QtWidgets.QVBoxLayout(self.layoutWidget)
        self.verticalLayout_8.setSpacing(15)
        self.verticalLayout_8.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout_8.setObjectName("verticalLayout_8")
        self.groupBox_18 = QtWidgets.QGroupBox(self.layoutWidget)
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        self.groupBox_18.setFont(font)
        self.groupBox_18.setStyleSheet("QGroupBox {\n"
"  \n"
"  border: 1px solid rgb(220, 220, 220);\n"
"  border-radius: 14px;\n"
"  padding: 4px;\n"
"  margin-top: 16px;\n"
"}\n"
"\n"
"QGroupBox::title {\n"
"  subcontrol-origin: margin;\n"
"  subcontrol-position: top left;\n"
"  left: 10px;\n"
"  padding-left: 3px;\n"
"  padding-right: 5px;\n"
"  padding-top: 8px;\n"
"  padding-bottom: 16px;\n"
"}\n"
"\n"
"QGroupBox::indicator {\n"
"  margin-left: 2px;\n"
"  height: 16px;\n"
"  width: 16px;\n"
"}")
        self.groupBox_18.setObjectName("groupBox_18")
        self.formLayoutWidget_9 = QtWidgets.QWidget(self.groupBox_18)
        self.formLayoutWidget_9.setGeometry(QtCore.QRect(10, 29, 401, 81))
        self.formLayoutWidget_9.setObjectName("formLayoutWidget_9")
        self.gridLayout_7 = QtWidgets.QGridLayout(self.formLayoutWidget_9)
        self.gridLayout_7.setContentsMargins(3, 4, 3, 4)
        self.gridLayout_7.setObjectName("gridLayout_7")
        self.label_39 = QtWidgets.QLabel(self.formLayoutWidget_9)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_39.sizePolicy().hasHeightForWidth())
        self.label_39.setSizePolicy(sizePolicy)
        self.label_39.setMinimumSize(QtCore.QSize(80, 0))
        self.label_39.setMaximumSize(QtCore.QSize(80, 16777215))
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(10)
        self.label_39.setFont(font)
        self.label_39.setObjectName("label_39")
        self.gridLayout_7.addWidget(self.label_39, 1, 0, 1, 1)
        self.platform_comboBox = QtWidgets.QComboBox(self.formLayoutWidget_9)
        self.platform_comboBox.setMinimumSize(QtCore.QSize(110, 33))
        self.platform_comboBox.setMaximumSize(QtCore.QSize(180, 33))
        self.platform_comboBox.setStyleSheet("QComboBox {\n"
"    background-color: rgb(27, 29, 35);\n"
"    border-radius: 5px;\n"
"    border: 2px solid rgb(27, 29, 35);\n"
"    padding: 1px 18px 1px 3px;\n"
"    padding-left: 10px;\n"
"    min-width: 6em;\n"
"}\n"
"QComboBox:hover{\n"
"    border: 2px solid rgb(64, 71, 88);\n"
"}\n"
"\n"
"QComboBox QAbstractItemView{\n"
"    \n"
"    selection-background-color: rgb(200, 5, 100);\n"
"    \n"
"\n"
"}\n"
"\n"
"QComboBox:editable {\n"
"\n"
"    background: rgb(0, 180, 95);\n"
"}\n"
"\n"
"QComboBox:on { /* shift the text when the popup opens */\n"
"    padding-top: 3px;\n"
"    padding-left: 4px;\n"
"    text-align: center;\n"
"}\n"
"\n"
"\n"
"\n"
"QComboBox:on { /* shift the text when the popup opens */\n"
"    padding-top: 2px;\n"
"    padding-left: 3px;\n"
"}\n"
"\n"
"QComboBox::drop-down {\n"
"    subcontrol-origin: padding;\n"
"    subcontrol-position: top right;\n"
"    width: 20px;\n"
"\n"
"    border-left-width: 2px;\n"
"    border-left-color: darkgray;\n"
"    border-left-style: solid; /* just a single line */\n"
"    border-top-right-radius: 3px; /* same radius as the QComboBox */\n"
"    border-bottom-right-radius: 3px;\n"
"}\n"
"\n"
"QComboBox::down-arrow {\n"
"image: url(:/16x16/icons/16x16/cil-arrow-bottom.png);\n"
"}\n"
"\n"
"QComboBox::down-arrow:on { /* shift the arrow when popup is open */\n"
"    top: 1px;\n"
"    left: 1px;\n"
"}\n"
"\n"
"")
        self.platform_comboBox.setObjectName("platform_comboBox")
        self.platform_comboBox.addItem("")
        self.platform_comboBox.addItem("")
        self.platform_comboBox.addItem("")
        self.platform_comboBox.addItem("")
        self.gridLayout_7.addWidget(self.platform_comboBox, 0, 1, 1, 1)
        self.Device_Number_comboBox = QtWidgets.QComboBox(self.formLayoutWidget_9)
        self.Device_Number_comboBox.setMinimumSize(QtCore.QSize(110, 33))
        self.Device_Number_comboBox.setMaximumSize(QtCore.QSize(180, 33))
        self.Device_Number_comboBox.setStyleSheet("QComboBox {\n"
"    background-color: rgb(27, 29, 35);\n"
"    border-radius: 5px;\n"
"    border: 2px solid rgb(27, 29, 35);\n"
"    padding: 1px 18px 1px 3px;\n"
"    padding-left: 10px;\n"
"    min-width: 6em;\n"
"}\n"
"QComboBox:hover{\n"
"    border: 2px solid rgb(64, 71, 88);\n"
"}\n"
"\n"
"QComboBox:disabled {\n"
"background-color: rgb(93, 83, 91);\n"
"}\n"
"\n"
"QComboBox QAbstractItemView{\n"
"    \n"
"    selection-background-color: rgb(200, 5, 100);\n"
"    \n"
"\n"
"}\n"
"\n"
"QComboBox:editable {\n"
"\n"
"    background: rgb(0, 180, 95);\n"
"}\n"
"\n"
"QComboBox:on { /* shift the text when the popup opens */\n"
"    padding-top: 3px;\n"
"    padding-left: 4px;\n"
"    text-align: center;\n"
"}\n"
"\n"
"\n"
"\n"
"QComboBox:on { /* shift the text when the popup opens */\n"
"    padding-top: 2px;\n"
"    padding-left: 3px;\n"
"}\n"
"\n"
"QComboBox::drop-down {\n"
"    subcontrol-origin: padding;\n"
"    subcontrol-position: top right;\n"
"    width: 20px;\n"
"\n"
"    border-left-width: 2px;\n"
"    border-left-color: darkgray;\n"
"    border-left-style: solid; /* just a single line */\n"
"    border-top-right-radius: 3px; /* same radius as the QComboBox */\n"
"    border-bottom-right-radius: 3px;\n"
"}\n"
"\n"
"QComboBox::down-arrow {\n"
"image: url(:/16x16/icons/16x16/cil-arrow-bottom.png);\n"
"}\n"
"\n"
"QComboBox::down-arrow:on { /* shift the arrow when popup is open */\n"
"    top: 1px;\n"
"    left: 1px;\n"
"}\n"
"\n"
"")
        self.Device_Number_comboBox.setObjectName("Device_Number_comboBox")
        self.Device_Number_comboBox.addItem("")
        self.gridLayout_7.addWidget(self.Device_Number_comboBox, 1, 1, 1, 1)
        self.label_38 = QtWidgets.QLabel(self.formLayoutWidget_9)
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(10)
        self.label_38.setFont(font)
        self.label_38.setObjectName("label_38")
        self.gridLayout_7.addWidget(self.label_38, 0, 0, 1, 1)
        self.Device_ID_checkBox = QtWidgets.QCheckBox(self.formLayoutWidget_9)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.Device_ID_checkBox.sizePolicy().hasHeightForWidth())
        self.Device_ID_checkBox.setSizePolicy(sizePolicy)
        self.Device_ID_checkBox.setMinimumSize(QtCore.QSize(120, 33))
        self.Device_ID_checkBox.setMaximumSize(QtCore.QSize(120, 33))
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(10)
        self.Device_ID_checkBox.setFont(font)
        self.Device_ID_checkBox.setStyleSheet("QCheckBox {\n"
"    spacing: 5px;\n"
"    margin-left: 10px;\n"
"}\n"
"\n"
"QCheckBox::indicator {\n"
"    width: 15px;\n"
"    height: 15px;\n"
"}\n"
"\n"
"QCheckBox::indicator:unchecked {\n"
"border-image: url(:/16x16/icons/16x16/cil-circle.png);\n"
"\n"
"}\n"
"\n"
"QCheckBox:disabled {\n"
"background-color: rgb(93, 83, 91);\n"
"}\n"
"\n"
"QCheckBox::indicator:checked {\n"
"    border-image: url(:/20x20/icons/20x20/cil-check.png);\n"
"}\n"
"\n"
"\n"
"")
        self.Device_ID_checkBox.setShortcut("")
        self.Device_ID_checkBox.setObjectName("Device_ID_checkBox")
        self.gridLayout_7.addWidget(self.Device_ID_checkBox, 1, 2, 1, 1)
        self.verticalLayout_8.addWidget(self.groupBox_18)
        self.groupBox_4 = QtWidgets.QGroupBox(self.layoutWidget)
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        self.groupBox_4.setFont(font)
        self.groupBox_4.setStyleSheet("QGroupBox {\n"
"  \n"
"  border: 1px solid rgb(220, 220, 220);\n"
"  border-radius: 14px;\n"
"  padding: 4px;\n"
"  margin-top: 16px;\n"
"}\n"
"\n"
"QGroupBox::title {\n"
"  subcontrol-origin: margin;\n"
"  subcontrol-position: top left;\n"
"  left: 10px;\n"
"  padding-left: 3px;\n"
"  padding-right: 5px;\n"
"  padding-top: 8px;\n"
"  padding-bottom: 16px;\n"
"}\n"
"\n"
"QGroupBox::indicator {\n"
"  margin-left: 2px;\n"
"  height: 16px;\n"
"  width: 16px;\n"
"}\n"
"\n"
"")
        self.groupBox_4.setObjectName("groupBox_4")
        self.formLayoutWidget = QtWidgets.QWidget(self.groupBox_4)
        self.formLayoutWidget.setGeometry(QtCore.QRect(9, 30, 561, 81))
        self.formLayoutWidget.setObjectName("formLayoutWidget")
        self.gridLayout_8 = QtWidgets.QGridLayout(self.formLayoutWidget)
        self.gridLayout_8.setContentsMargins(3, 4, 3, 4)
        self.gridLayout_8.setObjectName("gridLayout_8")
        self.protein_forcefield_comboBox = QtWidgets.QComboBox(self.formLayoutWidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.protein_forcefield_comboBox.sizePolicy().hasHeightForWidth())
        self.protein_forcefield_comboBox.setSizePolicy(sizePolicy)
        self.protein_forcefield_comboBox.setMinimumSize(QtCore.QSize(110, 33))
        self.protein_forcefield_comboBox.setMaximumSize(QtCore.QSize(180, 33))
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(10)
        self.protein_forcefield_comboBox.setFont(font)
        self.protein_forcefield_comboBox.setStyleSheet("QComboBox {\n"
"    background-color: rgb(27, 29, 35);\n"
"    border-radius: 5px;\n"
"    border: 2px solid rgb(27, 29, 35);\n"
"    padding: 1px 18px 1px 3px;\n"
"    padding-left: 10px;\n"
"    min-width: 6em;\n"
"}\n"
"QComboBox:hover{\n"
"    border: 2px solid rgb(64, 71, 88);\n"
"}\n"
"\n"
"QComboBox QAbstractItemView{\n"
"    \n"
"    selection-background-color: rgb(200, 5, 100);\n"
"    \n"
"\n"
"}\n"
"\n"
"QComboBox:editable {\n"
"\n"
"    background: rgb(0, 180, 95);\n"
"}\n"
"\n"
"QComboBox:on { /* shift the text when the popup opens */\n"
"    padding-top: 3px;\n"
"    padding-left: 4px;\n"
"    text-align: center;\n"
"}\n"
"\n"
"\n"
"\n"
"QComboBox:on { /* shift the text when the popup opens */\n"
"    padding-top: 2px;\n"
"    padding-left: 3px;\n"
"}\n"
"\n"
"QComboBox::drop-down {\n"
"    subcontrol-origin: padding;\n"
"    subcontrol-position: top right;\n"
"    width: 20px;\n"
"\n"
"    border-left-width: 2px;\n"
"    border-left-color: darkgray;\n"
"    border-left-style: solid; /* just a single line */\n"
"    border-top-right-radius: 3px; /* same radius as the QComboBox */\n"
"    border-bottom-right-radius: 3px;\n"
"}\n"
"\n"
"QComboBox::down-arrow {\n"
"image: url(:/16x16/icons/16x16/cil-arrow-bottom.png);\n"
"}\n"
"\n"
"QComboBox::down-arrow:on { /* shift the arrow when popup is open */\n"
"    top: 1px;\n"
"    left: 1px;\n"
"}\n"
"\n"
"")
        self.protein_forcefield_comboBox.setObjectName("protein_forcefield_comboBox")
        self.protein_forcefield_comboBox.addItem("")
        self.protein_forcefield_comboBox.addItem("")
        self.protein_forcefield_comboBox.addItem("")
        self.protein_forcefield_comboBox.addItem("")
        self.protein_forcefield_comboBox.addItem("")
        self.protein_forcefield_comboBox.addItem("")
        self.gridLayout_8.addWidget(self.protein_forcefield_comboBox, 0, 1, 1, 1)
        self.label_10 = QtWidgets.QLabel(self.formLayoutWidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_10.sizePolicy().hasHeightForWidth())
        self.label_10.setSizePolicy(sizePolicy)
        self.label_10.setMinimumSize(QtCore.QSize(80, 0))
        self.label_10.setMaximumSize(QtCore.QSize(80, 16777215))
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(10)
        self.label_10.setFont(font)
        self.label_10.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignVCenter)
        self.label_10.setObjectName("label_10")
        self.gridLayout_8.addWidget(self.label_10, 0, 0, 1, 1)
        self.label_11 = QtWidgets.QLabel(self.formLayoutWidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_11.sizePolicy().hasHeightForWidth())
        self.label_11.setSizePolicy(sizePolicy)
        self.label_11.setMinimumSize(QtCore.QSize(80, 0))
        self.label_11.setMaximumSize(QtCore.QSize(80, 16777215))
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(10)
        self.label_11.setFont(font)
        self.label_11.setObjectName("label_11")
        self.gridLayout_8.addWidget(self.label_11, 1, 0, 1, 1)
        self.water_forcefield_comboBox = QtWidgets.QComboBox(self.formLayoutWidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.water_forcefield_comboBox.sizePolicy().hasHeightForWidth())
        self.water_forcefield_comboBox.setSizePolicy(sizePolicy)
        self.water_forcefield_comboBox.setMinimumSize(QtCore.QSize(110, 33))
        self.water_forcefield_comboBox.setMaximumSize(QtCore.QSize(180, 33))
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(10)
        self.water_forcefield_comboBox.setFont(font)
        self.water_forcefield_comboBox.setStyleSheet("QComboBox {\n"
"    background-color: rgb(27, 29, 35);\n"
"    border-radius: 5px;\n"
"    border: 2px solid rgb(27, 29, 35);\n"
"    padding: 1px 18px 1px 3px;\n"
"    padding-left: 10px;\n"
"    min-width: 6em;\n"
"}\n"
"QComboBox:hover{\n"
"    border: 2px solid rgb(64, 71, 88);\n"
"}\n"
"\n"
"QComboBox QAbstractItemView{\n"
"    \n"
"    selection-background-color: rgb(200, 5, 100);\n"
"    \n"
"\n"
"}\n"
"\n"
"QComboBox:editable {\n"
"\n"
"    background: rgb(0, 180, 95);\n"
"}\n"
"\n"
"QComboBox:on { /* shift the text when the popup opens */\n"
"    padding-top: 3px;\n"
"    padding-left: 4px;\n"
"    text-align: center;\n"
"}\n"
"\n"
"\n"
"\n"
"QComboBox:on { /* shift the text when the popup opens */\n"
"    padding-top: 2px;\n"
"    padding-left: 3px;\n"
"}\n"
"\n"
"QComboBox::drop-down {\n"
"    subcontrol-origin: padding;\n"
"    subcontrol-position: top right;\n"
"    width: 20px;\n"
"\n"
"    border-left-width: 2px;\n"
"    border-left-color: darkgray;\n"
"    border-left-style: solid; /* just a single line */\n"
"    border-top-right-radius: 3px; /* same radius as the QComboBox */\n"
"    border-bottom-right-radius: 3px;\n"
"}\n"
"\n"
"QComboBox::down-arrow {\n"
"image: url(:/16x16/icons/16x16/cil-arrow-bottom.png);\n"
"}\n"
"\n"
"QComboBox::down-arrow:on { /* shift the arrow when popup is open */\n"
"    top: 1px;\n"
"    left: 1px;\n"
"}\n"
"\n"
"")
        self.water_forcefield_comboBox.setObjectName("water_forcefield_comboBox")
        self.water_forcefield_comboBox.addItem("")
        self.water_forcefield_comboBox.addItem("")
        self.water_forcefield_comboBox.addItem("")
        self.water_forcefield_comboBox.addItem("")
        self.gridLayout_8.addWidget(self.water_forcefield_comboBox, 1, 1, 1, 1)
        self.label_42 = QtWidgets.QLabel(self.formLayoutWidget)
        self.label_42.setMinimumSize(QtCore.QSize(80, 33))
        self.label_42.setMaximumSize(QtCore.QSize(80, 33))
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(10)
        self.label_42.setFont(font)
        self.label_42.setStyleSheet("QLabel {\n"
"    margin-left: 10px;\n"
"}")
        self.label_42.setObjectName("label_42")
        self.gridLayout_8.addWidget(self.label_42, 1, 2, 1, 1)
        self.water_padding_lineEdit = QtWidgets.QLineEdit(self.formLayoutWidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.water_padding_lineEdit.sizePolicy().hasHeightForWidth())
        self.water_padding_lineEdit.setSizePolicy(sizePolicy)
        self.water_padding_lineEdit.setMinimumSize(QtCore.QSize(150, 33))
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(10)
        self.water_padding_lineEdit.setFont(font)
        self.water_padding_lineEdit.setStyleSheet("QLineEdit {\n"
"    background-color: rgb(27, 29, 35);\n"
"    border-radius: 5px;\n"
"    padding-left: 10px;\n"
"\n"
"}\n"
"QLineEdit:hover {\n"
"    border: 2px solid rgb(64, 71, 88);\n"
"}\n"
"QLineEdit:focus {\n"
"    border: 2px solid rgb(91, 101, 124);\n"
"}")
        self.water_padding_lineEdit.setObjectName("water_padding_lineEdit")
        self.gridLayout_8.addWidget(self.water_padding_lineEdit, 1, 3, 1, 1)
        self.verticalLayout_8.addWidget(self.groupBox_4)
        self.tabWidget.addTab(self.tab, "")
        self.tab_2 = QtWidgets.QWidget()
        self.tab_2.setObjectName("tab_2")
        self.verticalLayout_20 = QtWidgets.QVBoxLayout(self.tab_2)
        self.verticalLayout_20.setObjectName("verticalLayout_20")
        self.groupBox_5 = QtWidgets.QGroupBox(self.tab_2)
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(10)
        self.groupBox_5.setFont(font)
        self.groupBox_5.setStyleSheet("QGroupBox {\n"
"  \n"
"  border: 1px solid rgb(220, 220, 220);\n"
"  border-radius: 10px;\n"
"  padding: 4px;\n"
"  margin-top: 16px;\n"
"}\n"
"\n"
"QGroupBox::title {\n"
"  subcontrol-origin: margin;\n"
"  subcontrol-position: top left;\n"
"  left: 10px;\n"
"  padding-left: 3px;\n"
"  padding-right: 5px;\n"
"  padding-top: 8px;\n"
"  padding-bottom: 16px;\n"
"}\n"
"\n"
"QGroupBox::indicator {\n"
"  margin-left: 2px;\n"
"  height: 16px;\n"
"  width: 16px;\n"
"}")
        self.groupBox_5.setObjectName("groupBox_5")
        self.verticalLayout_18 = QtWidgets.QVBoxLayout(self.groupBox_5)
        self.verticalLayout_18.setObjectName("verticalLayout_18")
        self.groupBox_3 = QtWidgets.QGroupBox(self.groupBox_5)
        self.groupBox_3.setMinimumSize(QtCore.QSize(0, 150))
        self.groupBox_3.setMaximumSize(QtCore.QSize(16777215, 150))
        self.groupBox_3.setObjectName("groupBox_3")
        self.horizontalLayout_19 = QtWidgets.QHBoxLayout(self.groupBox_3)
        self.horizontalLayout_19.setObjectName("horizontalLayout_19")
        self.gridLayout = QtWidgets.QGridLayout()
        self.gridLayout.setContentsMargins(3, 4, 3, 4)
        self.gridLayout.setHorizontalSpacing(6)
        self.gridLayout.setVerticalSpacing(0)
        self.gridLayout.setObjectName("gridLayout")
        self.integrator_kind_comboBox = QtWidgets.QComboBox(self.groupBox_3)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.integrator_kind_comboBox.sizePolicy().hasHeightForWidth())
        self.integrator_kind_comboBox.setSizePolicy(sizePolicy)
        self.integrator_kind_comboBox.setMinimumSize(QtCore.QSize(110, 38))
        self.integrator_kind_comboBox.setMaximumSize(QtCore.QSize(240, 38))
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(10)
        self.integrator_kind_comboBox.setFont(font)
        self.integrator_kind_comboBox.setStyleSheet("QComboBox {\n"
"    background-color: rgb(27, 29, 35);\n"
"    border-radius: 5px;\n"
"    border: 2px solid rgb(27, 29, 35);\n"
"    padding: 1px 18px 1px 3px;\n"
"    padding-left: 10px;\n"
"    min-width: 6em;\n"
"}\n"
"QComboBox:hover{\n"
"    border: 2px solid rgb(64, 71, 88);\n"
"}\n"
"\n"
"QComboBox QAbstractItemView{\n"
"    \n"
"    selection-background-color: rgb(200, 5, 100);\n"
"    \n"
"\n"
"}\n"
"\n"
"QComboBox:editable {\n"
"\n"
"    background: rgb(0, 180, 95);\n"
"}\n"
"\n"
"QComboBox:on { /* shift the text when the popup opens */\n"
"    padding-top: 3px;\n"
"    padding-left: 4px;\n"
"    text-align: center;\n"
"}\n"
"\n"
"\n"
"\n"
"QComboBox:on { /* shift the text when the popup opens */\n"
"    padding-top: 2px;\n"
"    padding-left: 3px;\n"
"}\n"
"\n"
"QComboBox::drop-down {\n"
"    subcontrol-origin: padding;\n"
"    subcontrol-position: top right;\n"
"    width: 20px;\n"
"\n"
"    border-left-width: 2px;\n"
"    border-left-color: darkgray;\n"
"    border-left-style: solid; /* just a single line */\n"
"    border-top-right-radius: 3px; /* same radius as the QComboBox */\n"
"    border-bottom-right-radius: 3px;\n"
"}\n"
"\n"
"QComboBox::down-arrow {\n"
"image: url(:/16x16/icons/16x16/cil-arrow-bottom.png);\n"
"}\n"
"\n"
"QComboBox::down-arrow:on { /* shift the arrow when popup is open */\n"
"    top: 1px;\n"
"    left: 1px;\n"
"}\n"
"\n"
"")
        self.integrator_kind_comboBox.setObjectName("integrator_kind_comboBox")
        self.integrator_kind_comboBox.addItem("")
        self.integrator_kind_comboBox.addItem("")
        self.integrator_kind_comboBox.addItem("")
        self.gridLayout.addWidget(self.integrator_kind_comboBox, 1, 1, 2, 1)
        self.integrator_time_step = QtWidgets.QTextEdit(self.groupBox_3)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.integrator_time_step.sizePolicy().hasHeightForWidth())
        self.integrator_time_step.setSizePolicy(sizePolicy)
        self.integrator_time_step.setMinimumSize(QtCore.QSize(0, 38))
        self.integrator_time_step.setMaximumSize(QtCore.QSize(240, 38))
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(10)
        self.integrator_time_step.setFont(font)
        self.integrator_time_step.setStyleSheet("QTextEdit {\n"
"    background-color: rgb(27, 29, 35);\n"
"    border-radius: 5px;\n"
"    padding: 10px;\n"
"    \n"
"    color: rgb(255, 255, 255);\n"
"}\n"
"QTextEdit:hover {\n"
"    border: 2px solid rgb(64, 71, 88);\n"
"}\n"
"QTextEdit:focus {\n"
"    border: 2px solid rgb(91, 101, 124);\n"
"}")
        self.integrator_time_step.setObjectName("integrator_time_step")
        self.gridLayout.addWidget(self.integrator_time_step, 3, 1, 1, 1)
        self.label_12 = QtWidgets.QLabel(self.groupBox_3)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_12.sizePolicy().hasHeightForWidth())
        self.label_12.setSizePolicy(sizePolicy)
        self.label_12.setMinimumSize(QtCore.QSize(85, 0))
        self.label_12.setMaximumSize(QtCore.QSize(85, 16777215))
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(10)
        self.label_12.setFont(font)
        self.label_12.setObjectName("label_12")
        self.gridLayout.addWidget(self.label_12, 3, 0, 1, 1)
        self.label_9 = QtWidgets.QLabel(self.groupBox_3)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_9.sizePolicy().hasHeightForWidth())
        self.label_9.setSizePolicy(sizePolicy)
        self.label_9.setMinimumSize(QtCore.QSize(85, 0))
        self.label_9.setMaximumSize(QtCore.QSize(85, 16777215))
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(10)
        self.label_9.setFont(font)
        self.label_9.setObjectName("label_9")
        self.gridLayout.addWidget(self.label_9, 1, 0, 1, 1)
        self.Additional_Integrators_checkBox = QtWidgets.QCheckBox(self.groupBox_3)
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(10)
        self.Additional_Integrators_checkBox.setFont(font)
        self.Additional_Integrators_checkBox.setStyleSheet("QCheckBox {\n"
"    spacing: 5px;\n"
"    margin-left: 15px;\n"
"}\n"
"\n"
"QCheckBox::indicator {\n"
"    width: 15px;\n"
"    height: 15px;\n"
"}\n"
"\n"
"QCheckBox::indicator:unchecked {\n"
"border-image: url(:/16x16/icons/16x16/cil-circle.png);\n"
"\n"
"}\n"
"\n"
"\n"
"QCheckBox::indicator:checked {\n"
"    border-image: url(:/20x20/icons/20x20/cil-check.png);\n"
"}\n"
"\n"
"\n"
"")
        self.Additional_Integrators_checkBox.setChecked(True)
        self.Additional_Integrators_checkBox.setObjectName("Additional_Integrators_checkBox")
        self.gridLayout.addWidget(self.Additional_Integrators_checkBox, 1, 2, 1, 1)
        self.integrator_time_step_unit = QtWidgets.QComboBox(self.groupBox_3)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.integrator_time_step_unit.sizePolicy().hasHeightForWidth())
        self.integrator_time_step_unit.setSizePolicy(sizePolicy)
        self.integrator_time_step_unit.setMaximumSize(QtCore.QSize(16777215, 38))
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(10)
        self.integrator_time_step_unit.setFont(font)
        self.integrator_time_step_unit.setStyleSheet("QComboBox {\n"
"    background-color: rgb(27, 29, 35);\n"
"    border-radius: 5px;\n"
"    border: 2px solid rgb(27, 29, 35);\n"
"    padding: 1px 18px 1px 3px;\n"
"    padding-left: 10px;\n"
"    min-width: 6em;\n"
"}\n"
"QComboBox:hover{\n"
"    border: 2px solid rgb(64, 71, 88);\n"
"}\n"
"\n"
"QComboBox QAbstractItemView{\n"
"    \n"
"    selection-background-color: rgb(200, 5, 100);\n"
"    \n"
"\n"
"}\n"
"\n"
"QComboBox:editable {\n"
"\n"
"    background: rgb(0, 180, 95);\n"
"}\n"
"\n"
"QComboBox:on { /* shift the text when the popup opens */\n"
"    padding-top: 3px;\n"
"    padding-left: 4px;\n"
"    text-align: center;\n"
"}\n"
"\n"
"\n"
"\n"
"QComboBox:on { /* shift the text when the popup opens */\n"
"    padding-top: 2px;\n"
"    padding-left: 3px;\n"
"}\n"
"\n"
"QComboBox::drop-down {\n"
"    subcontrol-origin: padding;\n"
"    subcontrol-position: top right;\n"
"    width: 20px;\n"
"\n"
"    border-left-width: 2px;\n"
"    border-left-color: darkgray;\n"
"    border-left-style: solid; /* just a single line */\n"
"    border-top-right-radius: 3px; /* same radius as the QComboBox */\n"
"    border-bottom-right-radius: 3px;\n"
"}\n"
"\n"
"QComboBox::down-arrow {\n"
"image: url(:/16x16/icons/16x16/cil-arrow-bottom.png);\n"
"}\n"
"\n"
"QComboBox::down-arrow:on { /* shift the arrow when popup is open */\n"
"    top: 1px;\n"
"    left: 1px;\n"
"}\n"
"\n"
"")
        self.integrator_time_step_unit.setObjectName("integrator_time_step_unit")
        self.integrator_time_step_unit.addItem("")
        self.gridLayout.addWidget(self.integrator_time_step_unit, 3, 2, 1, 1)
        self.horizontalLayout_19.addLayout(self.gridLayout)
        self.verticalLayout_18.addWidget(self.groupBox_3)
        self.Additional_Integrator_groupBox = QtWidgets.QGroupBox(self.groupBox_5)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.Additional_Integrator_groupBox.sizePolicy().hasHeightForWidth())
        self.Additional_Integrator_groupBox.setSizePolicy(sizePolicy)
        self.Additional_Integrator_groupBox.setMinimumSize(QtCore.QSize(0, 150))
        self.Additional_Integrator_groupBox.setMaximumSize(QtCore.QSize(16777215, 150))
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(10)
        self.Additional_Integrator_groupBox.setFont(font)
        self.Additional_Integrator_groupBox.setStyleSheet("QGroupBox {\n"
"  \n"
"  border: 1px solid rgb(220, 220, 220);\n"
"  border-radius: 10px;\n"
"  padding: 4px;\n"
"  margin-top: 16px;\n"
"}\n"
"\n"
"QGroupBox::title {\n"
"  subcontrol-origin: margin;\n"
"  subcontrol-position: top left;\n"
"  left: 10px;\n"
"  padding-left: 3px;\n"
"  padding-right: 5px;\n"
"  padding-top: 8px;\n"
"  padding-bottom: 16px;\n"
"}\n"
"\n"
"QGroupBox::indicator {\n"
"  margin-left: 2px;\n"
"  height: 16px;\n"
"  width: 16px;\n"
"}\n"
"\n"
"QWidget:disabled {\n"
"  background-color: #19232D;\n"
"  color: #787878;\n"
"  selection-background-color: #14506E;\n"
"  selection-color: #787878;\n"
"}\n"
"\n"
"")
        self.Additional_Integrator_groupBox.setObjectName("Additional_Integrator_groupBox")
        self.formLayoutWidget_3 = QtWidgets.QWidget(self.Additional_Integrator_groupBox)
        self.formLayoutWidget_3.setGeometry(QtCore.QRect(9, 30, 341, 90))
        self.formLayoutWidget_3.setObjectName("formLayoutWidget_3")
        self.formLayout_3 = QtWidgets.QFormLayout(self.formLayoutWidget_3)
        self.formLayout_3.setContentsMargins(3, 4, 3, 4)
        self.formLayout_3.setVerticalSpacing(6)
        self.formLayout_3.setObjectName("formLayout_3")
        self.label_13 = QtWidgets.QLabel(self.formLayoutWidget_3)
        self.label_13.setMinimumSize(QtCore.QSize(85, 0))
        self.label_13.setMaximumSize(QtCore.QSize(85, 16777215))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.label_13.setFont(font)
        self.label_13.setObjectName("label_13")
        self.formLayout_3.setWidget(0, QtWidgets.QFormLayout.LabelRole, self.label_13)
        self.label_14 = QtWidgets.QLabel(self.formLayoutWidget_3)
        self.label_14.setMinimumSize(QtCore.QSize(85, 0))
        self.label_14.setMaximumSize(QtCore.QSize(85, 16777215))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.label_14.setFont(font)
        self.label_14.setObjectName("label_14")
        self.formLayout_3.setWidget(1, QtWidgets.QFormLayout.LabelRole, self.label_14)
        self.friction_textEdit = QtWidgets.QTextEdit(self.formLayoutWidget_3)
        self.friction_textEdit.setMinimumSize(QtCore.QSize(0, 38))
        self.friction_textEdit.setMaximumSize(QtCore.QSize(240, 38))
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(10)
        self.friction_textEdit.setFont(font)
        self.friction_textEdit.setStyleSheet("QTextEdit {\n"
"    background-color: rgb(27, 29, 35);\n"
"    border-radius: 5px;\n"
"    padding: 10px;\n"
"    \n"
"    color: rgb(255, 255, 255);\n"
"}\n"
"QTextEdit:hover {\n"
"    border: 2px solid rgb(64, 71, 88);\n"
"}\n"
"QTextEdit:focus {\n"
"    border: 2px solid rgb(91, 101, 124);\n"
"}")
        self.friction_textEdit.setObjectName("friction_textEdit")
        self.formLayout_3.setWidget(0, QtWidgets.QFormLayout.FieldRole, self.friction_textEdit)
        self.temperature_textEdit = QtWidgets.QTextEdit(self.formLayoutWidget_3)
        self.temperature_textEdit.setMinimumSize(QtCore.QSize(0, 38))
        self.temperature_textEdit.setMaximumSize(QtCore.QSize(240, 38))
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(10)
        self.temperature_textEdit.setFont(font)
        self.temperature_textEdit.setStyleSheet("QTextEdit {\n"
"    background-color: rgb(27, 29, 35);\n"
"    border-radius: 5px;\n"
"    padding: 10px;\n"
"    \n"
"    color: rgb(255, 255, 255);\n"
"}\n"
"QTextEdit:hover {\n"
"    border: 2px solid rgb(64, 71, 88);\n"
"}\n"
"QTextEdit:focus {\n"
"    border: 2px solid rgb(91, 101, 124);\n"
"}")
        self.temperature_textEdit.setObjectName("temperature_textEdit")
        self.formLayout_3.setWidget(1, QtWidgets.QFormLayout.FieldRole, self.temperature_textEdit)
        self.verticalLayout_18.addWidget(self.Additional_Integrator_groupBox)
        spacerItem1 = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.verticalLayout_18.addItem(spacerItem1)
        self.verticalLayout_20.addWidget(self.groupBox_5)
        self.tabWidget.addTab(self.tab_2, "")
        self.tab_3 = QtWidgets.QWidget()
        self.tab_3.setObjectName("tab_3")
        self.groupBox_7 = QtWidgets.QGroupBox(self.tab_3)
        self.groupBox_7.setGeometry(QtCore.QRect(10, 10, 571, 251))
        font = QtGui.QFont()
        font.setFamily("Segoe UI Semibold")
        self.groupBox_7.setFont(font)
        self.groupBox_7.setStyleSheet("QGroupBox {\n"
"  \n"
"  border: 1px solid rgb(220, 220, 220);\n"
"  border-radius: 10px;\n"
"  padding: 4px;\n"
"  margin-top: 16px;\n"
"}\n"
"\n"
"QGroupBox::title {\n"
"  subcontrol-origin: margin;\n"
"  subcontrol-position: top left;\n"
"  left: 10px;\n"
"  padding-left: 3px;\n"
"  padding-right: 5px;\n"
"  padding-top: 8px;\n"
"  padding-bottom: 16px;\n"
"}\n"
"\n"
"QGroupBox::indicator {\n"
"  margin-left: 2px;\n"
"  height: 16px;\n"
"  width: 16px;\n"
"}")
        self.groupBox_7.setObjectName("groupBox_7")
        self.formLayoutWidget_2 = QtWidgets.QWidget(self.groupBox_7)
        self.formLayoutWidget_2.setGeometry(QtCore.QRect(10, 29, 551, 211))
        self.formLayoutWidget_2.setObjectName("formLayoutWidget_2")
        self.gridLayout_5 = QtWidgets.QGridLayout(self.formLayoutWidget_2)
        self.gridLayout_5.setContentsMargins(3, 4, 3, 4)
        self.gridLayout_5.setObjectName("gridLayout_5")
        self.system_constraints_comboBox = QtWidgets.QComboBox(self.formLayoutWidget_2)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.system_constraints_comboBox.sizePolicy().hasHeightForWidth())
        self.system_constraints_comboBox.setSizePolicy(sizePolicy)
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(10)
        self.system_constraints_comboBox.setFont(font)
        self.system_constraints_comboBox.setStyleSheet("QComboBox {\n"
"    background-color: rgb(27, 29, 35);\n"
"    border-radius: 5px;\n"
"    border: 2px solid rgb(27, 29, 35);\n"
"    padding: 1px 18px 1px 3px;\n"
"    padding-left: 10px;\n"
"    min-width: 6em;\n"
"}\n"
"QComboBox:hover{\n"
"    border: 2px solid rgb(64, 71, 88);\n"
"}\n"
"\n"
"QComboBox QAbstractItemView{\n"
"    \n"
"    selection-background-color: rgb(200, 5, 100);\n"
"    \n"
"\n"
"}\n"
"\n"
"QComboBox:editable {\n"
"\n"
"    background: rgb(0, 180, 95);\n"
"}\n"
"\n"
"QComboBox:on { /* shift the text when the popup opens */\n"
"    padding-top: 3px;\n"
"    padding-left: 4px;\n"
"    text-align: center;\n"
"}\n"
"\n"
"\n"
"\n"
"QComboBox:on { /* shift the text when the popup opens */\n"
"    padding-top: 2px;\n"
"    padding-left: 3px;\n"
"}\n"
"\n"
"QComboBox::drop-down {\n"
"    subcontrol-origin: padding;\n"
"    subcontrol-position: top right;\n"
"    width: 20px;\n"
"\n"
"    border-left-width: 2px;\n"
"    border-left-color: darkgray;\n"
"    border-left-style: solid; /* just a single line */\n"
"    border-top-right-radius: 3px; /* same radius as the QComboBox */\n"
"    border-bottom-right-radius: 3px;\n"
"}\n"
"\n"
"QComboBox::down-arrow {\n"
"image: url(:/16x16/icons/16x16/cil-arrow-bottom.png);\n"
"}\n"
"\n"
"QComboBox::down-arrow:on { /* shift the arrow when popup is open */\n"
"    top: 1px;\n"
"    left: 1px;\n"
"}\n"
"\n"
"")
        self.system_constraints_comboBox.setObjectName("system_constraints_comboBox")
        self.system_constraints_comboBox.addItem("")
        self.system_constraints_comboBox.addItem("")
        self.system_constraints_comboBox.addItem("")
        self.system_constraints_comboBox.addItem("")
        self.gridLayout_5.addWidget(self.system_constraints_comboBox, 1, 1, 1, 1)
        self.rigid_water_checkBox = QtWidgets.QCheckBox(self.formLayoutWidget_2)
        self.rigid_water_checkBox.setEnabled(True)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.rigid_water_checkBox.sizePolicy().hasHeightForWidth())
        self.rigid_water_checkBox.setSizePolicy(sizePolicy)
        self.rigid_water_checkBox.setStyleSheet("QCheckBox {\n"
"    spacing: 5px;\n"
"    margin-left: 15px;\n"
"}\n"
"\n"
"QCheckBox::indicator {\n"
"    width: 15px;\n"
"    height: 15px;\n"
"}\n"
"\n"
"QCheckBox::indicator:unchecked {\n"
"border-image: url(:/16x16/icons/16x16/cil-circle.png);\n"
"\n"
"}\n"
"\n"
"\n"
"QCheckBox::indicator:checked {\n"
"    border-image: url(:/20x20/icons/20x20/cil-check.png);\n"
"}\n"
"\n"
"\n"
"")
        self.rigid_water_checkBox.setText("")
        self.rigid_water_checkBox.setChecked(True)
        self.rigid_water_checkBox.setTristate(False)
        self.rigid_water_checkBox.setObjectName("rigid_water_checkBox")
        self.gridLayout_5.addWidget(self.rigid_water_checkBox, 2, 1, 1, 1)
        self.label_17 = QtWidgets.QLabel(self.formLayoutWidget_2)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_17.sizePolicy().hasHeightForWidth())
        self.label_17.setSizePolicy(sizePolicy)
        self.label_17.setMaximumSize(QtCore.QSize(160, 16777215))
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(10)
        self.label_17.setFont(font)
        self.label_17.setObjectName("label_17")
        self.gridLayout_5.addWidget(self.label_17, 2, 0, 1, 1)
        self.label_18 = QtWidgets.QLabel(self.formLayoutWidget_2)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_18.sizePolicy().hasHeightForWidth())
        self.label_18.setSizePolicy(sizePolicy)
        self.label_18.setMaximumSize(QtCore.QSize(160, 16777215))
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(10)
        self.label_18.setFont(font)
        self.label_18.setObjectName("label_18")
        self.gridLayout_5.addWidget(self.label_18, 3, 0, 1, 1)
        self.label_15 = QtWidgets.QLabel(self.formLayoutWidget_2)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_15.sizePolicy().hasHeightForWidth())
        self.label_15.setSizePolicy(sizePolicy)
        self.label_15.setMaximumSize(QtCore.QSize(160, 16777215))
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(10)
        self.label_15.setFont(font)
        self.label_15.setObjectName("label_15")
        self.gridLayout_5.addWidget(self.label_15, 0, 0, 1, 1)
        self.label_16 = QtWidgets.QLabel(self.formLayoutWidget_2)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_16.sizePolicy().hasHeightForWidth())
        self.label_16.setSizePolicy(sizePolicy)
        self.label_16.setMaximumSize(QtCore.QSize(160, 16777215))
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(10)
        self.label_16.setFont(font)
        self.label_16.setObjectName("label_16")
        self.gridLayout_5.addWidget(self.label_16, 1, 0, 1, 1)
        self.nonbounded_CutOff_textEdit = QtWidgets.QTextEdit(self.formLayoutWidget_2)
        self.nonbounded_CutOff_textEdit.setMinimumSize(QtCore.QSize(175, 38))
        self.nonbounded_CutOff_textEdit.setMaximumSize(QtCore.QSize(180, 38))
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(10)
        self.nonbounded_CutOff_textEdit.setFont(font)
        self.nonbounded_CutOff_textEdit.setStyleSheet("QTextEdit {\n"
"    background-color: rgb(27, 29, 35);\n"
"    border-radius: 5px;\n"
"    padding: 10px;\n"
"}\n"
"QTextEdit:hover {\n"
"    border: 2px solid rgb(64, 71, 88);\n"
"}\n"
"QTextEdit:focus {\n"
"    border: 2px solid rgb(91, 101, 124);\n"
"}")
        self.nonbounded_CutOff_textEdit.setObjectName("nonbounded_CutOff_textEdit")
        self.gridLayout_5.addWidget(self.nonbounded_CutOff_textEdit, 3, 1, 1, 1)
        self.nonBounded_Method_comboBox = QtWidgets.QComboBox(self.formLayoutWidget_2)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.nonBounded_Method_comboBox.sizePolicy().hasHeightForWidth())
        self.nonBounded_Method_comboBox.setSizePolicy(sizePolicy)
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(10)
        self.nonBounded_Method_comboBox.setFont(font)
        self.nonBounded_Method_comboBox.setStyleSheet("QComboBox {\n"
"    background-color: rgb(27, 29, 35);\n"
"    border-radius: 5px;\n"
"    border: 2px solid rgb(27, 29, 35);\n"
"    padding: 1px 18px 1px 3px;\n"
"    padding-left: 10px;\n"
"    min-width: 6em;\n"
"}\n"
"QComboBox:hover{\n"
"    border: 2px solid rgb(64, 71, 88);\n"
"}\n"
"\n"
"QComboBox QAbstractItemView{\n"
"    \n"
"    selection-background-color: rgb(200, 5, 100);\n"
"    \n"
"\n"
"}\n"
"\n"
"QComboBox:editable {\n"
"\n"
"    background: rgb(0, 180, 95);\n"
"}\n"
"\n"
"QComboBox:on { /* shift the text when the popup opens */\n"
"    padding-top: 3px;\n"
"    padding-left: 4px;\n"
"    text-align: center;\n"
"}\n"
"\n"
"\n"
"\n"
"QComboBox:on { /* shift the text when the popup opens */\n"
"    padding-top: 2px;\n"
"    padding-left: 3px;\n"
"}\n"
"\n"
"QComboBox::drop-down {\n"
"    subcontrol-origin: padding;\n"
"    subcontrol-position: top right;\n"
"    width: 20px;\n"
"\n"
"    border-left-width: 2px;\n"
"    border-left-color: darkgray;\n"
"    border-left-style: solid; /* just a single line */\n"
"    border-top-right-radius: 3px; /* same radius as the QComboBox */\n"
"    border-bottom-right-radius: 3px;\n"
"}\n"
"\n"
"QComboBox::down-arrow {\n"
"image: url(:/16x16/icons/16x16/cil-arrow-bottom.png);\n"
"}\n"
"\n"
"QComboBox::down-arrow:on { /* shift the arrow when popup is open */\n"
"    top: 1px;\n"
"    left: 1px;\n"
"}\n"
"\n"
"")
        self.nonBounded_Method_comboBox.setObjectName("nonBounded_Method_comboBox")
        self.nonBounded_Method_comboBox.addItem("")
        self.nonBounded_Method_comboBox.addItem("")
        self.nonBounded_Method_comboBox.addItem("")
        self.nonBounded_Method_comboBox.addItem("")
        self.nonBounded_Method_comboBox.addItem("")
        self.gridLayout_5.addWidget(self.nonBounded_Method_comboBox, 0, 1, 1, 1)
        self.switching_distance_textEdit = QtWidgets.QTextEdit(self.formLayoutWidget_2)
        self.switching_distance_textEdit.setMinimumSize(QtCore.QSize(0, 38))
        self.switching_distance_textEdit.setMaximumSize(QtCore.QSize(180, 38))
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(10)
        self.switching_distance_textEdit.setFont(font)
        self.switching_distance_textEdit.setStyleSheet("QTextEdit {\n"
"    background-color: rgb(27, 29, 35);\n"
"    border-radius: 5px;\n"
"    padding: 10px;\n"
"}\n"
"QTextEdit:hover {\n"
"    border: 2px solid rgb(64, 71, 88);\n"
"}\n"
"QTextEdit:focus {\n"
"    border: 2px solid rgb(91, 101, 124);\n"
"}")
        self.switching_distance_textEdit.setObjectName("switching_distance_textEdit")
        self.gridLayout_5.addWidget(self.switching_distance_textEdit, 4, 1, 1, 1)
        self.label_45 = QtWidgets.QLabel(self.formLayoutWidget_2)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_45.sizePolicy().hasHeightForWidth())
        self.label_45.setSizePolicy(sizePolicy)
        self.label_45.setMaximumSize(QtCore.QSize(160, 16777215))
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(10)
        self.label_45.setFont(font)
        self.label_45.setObjectName("label_45")
        self.gridLayout_5.addWidget(self.label_45, 4, 0, 1, 1)
        self.use_switching_checkBox = QtWidgets.QCheckBox(self.formLayoutWidget_2)
        self.use_switching_checkBox.setEnabled(True)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.use_switching_checkBox.sizePolicy().hasHeightForWidth())
        self.use_switching_checkBox.setSizePolicy(sizePolicy)
        self.use_switching_checkBox.setMinimumSize(QtCore.QSize(175, 0))
        self.use_switching_checkBox.setStyleSheet("QCheckBox {\n"
"    spacing: 5px;\n"
"    margin-left: 15px;\n"
"}\n"
"\n"
"QCheckBox::indicator {\n"
"    width: 15px;\n"
"    height: 15px;\n"
"}\n"
"\n"
"QCheckBox::indicator:unchecked {\n"
"border-image: url(:/16x16/icons/16x16/cil-circle.png);\n"
"\n"
"}\n"
"\n"
"\n"
"QCheckBox::indicator:checked {\n"
"    border-image: url(:/20x20/icons/20x20/cil-check.png);\n"
"}\n"
"\n"
"\n"
"")
        self.use_switching_checkBox.setChecked(True)
        self.use_switching_checkBox.setTristate(False)
        self.use_switching_checkBox.setObjectName("use_switching_checkBox")
        self.gridLayout_5.addWidget(self.use_switching_checkBox, 4, 2, 1, 1)
        self.tabWidget.addTab(self.tab_3, "")
        self.tab_4 = QtWidgets.QWidget()
        self.tab_4.setObjectName("tab_4")
        self.layoutWidget1 = QtWidgets.QWidget(self.tab_4)
        self.layoutWidget1.setGeometry(QtCore.QRect(10, 10, 581, 761))
        self.layoutWidget1.setObjectName("layoutWidget1")
        self.verticalLayout_11 = QtWidgets.QVBoxLayout(self.layoutWidget1)
        self.verticalLayout_11.setContentsMargins(0, -1, -1, 4)
        self.verticalLayout_11.setObjectName("verticalLayout_11")
        self.groupBox_8 = QtWidgets.QGroupBox(self.layoutWidget1)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.groupBox_8.sizePolicy().hasHeightForWidth())
        self.groupBox_8.setSizePolicy(sizePolicy)
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(10)
        self.groupBox_8.setFont(font)
        self.groupBox_8.setStyleSheet("QGroupBox {\n"
"  \n"
"  border: 1px solid rgb(220, 220, 220);\n"
"  border-radius: 10px;\n"
"  padding: 4px;\n"
"  margin-top: 16px;\n"
"}\n"
"\n"
"QGroupBox::title {\n"
"  subcontrol-origin: margin;\n"
"  subcontrol-position: top left;\n"
"  left: 10px;\n"
"  padding-left: 3px;\n"
"  padding-right: 5px;\n"
"  padding-top: 8px;\n"
"  padding-bottom: 16px;\n"
"}\n"
"\n"
"QGroupBox::indicator {\n"
"  margin-left: 2px;\n"
"  height: 16px;\n"
"  width: 16px;\n"
"}")
        self.groupBox_8.setObjectName("groupBox_8")
        self.groupBox_10 = QtWidgets.QGroupBox(self.groupBox_8)
        self.groupBox_10.setGeometry(QtCore.QRect(9, 100, 561, 121))
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(10)
        self.groupBox_10.setFont(font)
        self.groupBox_10.setObjectName("groupBox_10")
        self.formLayoutWidget_4 = QtWidgets.QWidget(self.groupBox_10)
        self.formLayoutWidget_4.setGeometry(QtCore.QRect(9, 30, 261, 81))
        self.formLayoutWidget_4.setObjectName("formLayoutWidget_4")
        self.formLayout_4 = QtWidgets.QFormLayout(self.formLayoutWidget_4)
        self.formLayout_4.setContentsMargins(3, 8, 3, 4)
        self.formLayout_4.setVerticalSpacing(9)
        self.formLayout_4.setObjectName("formLayout_4")
        self.label_20 = QtWidgets.QLabel(self.formLayoutWidget_4)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_20.sizePolicy().hasHeightForWidth())
        self.label_20.setSizePolicy(sizePolicy)
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(10)
        self.label_20.setFont(font)
        self.label_20.setObjectName("label_20")
        self.formLayout_4.setWidget(0, QtWidgets.QFormLayout.LabelRole, self.label_20)
        self.label_21 = QtWidgets.QLabel(self.formLayoutWidget_4)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_21.sizePolicy().hasHeightForWidth())
        self.label_21.setSizePolicy(sizePolicy)
        self.label_21.setMinimumSize(QtCore.QSize(107, 0))
        self.label_21.setMaximumSize(QtCore.QSize(107, 16777215))
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(10)
        self.label_21.setFont(font)
        self.label_21.setObjectName("label_21")
        self.formLayout_4.setWidget(1, QtWidgets.QFormLayout.LabelRole, self.label_21)
        self.minimize_checkBox = QtWidgets.QCheckBox(self.formLayoutWidget_4)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.minimize_checkBox.sizePolicy().hasHeightForWidth())
        self.minimize_checkBox.setSizePolicy(sizePolicy)
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(10)
        self.minimize_checkBox.setFont(font)
        self.minimize_checkBox.setStyleSheet("QCheckBox {\n"
"    spacing: 5px;\n"
"    margin-left: 15px;\n"
"}\n"
"\n"
"QCheckBox::indicator {\n"
"    width: 15px;\n"
"    height: 15px;\n"
"}\n"
"\n"
"QCheckBox::indicator:unchecked {\n"
"border-image: url(:/16x16/icons/16x16/cil-circle.png);\n"
"\n"
"}\n"
"\n"
"\n"
"QCheckBox::indicator:checked {\n"
"    border-image: url(:/20x20/icons/20x20/cil-check.png);\n"
"}\n"
"\n"
"\n"
"")
        self.minimize_checkBox.setText("")
        self.minimize_checkBox.setChecked(True)
        self.minimize_checkBox.setObjectName("minimize_checkBox")
        self.formLayout_4.setWidget(0, QtWidgets.QFormLayout.FieldRole, self.minimize_checkBox)
        self.Max_minimize_iter_textEdit = QtWidgets.QTextEdit(self.formLayoutWidget_4)
        self.Max_minimize_iter_textEdit.setMinimumSize(QtCore.QSize(0, 38))
        self.Max_minimize_iter_textEdit.setMaximumSize(QtCore.QSize(16777215, 38))
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(10)
        self.Max_minimize_iter_textEdit.setFont(font)
        self.Max_minimize_iter_textEdit.setStyleSheet("QTextEdit {\n"
"    background-color: rgb(27, 29, 35);\n"
"    border-radius: 5px;\n"
"    padding: 10px;\n"
"}\n"
"QTextEdit:hover {\n"
"    border: 2px solid rgb(64, 71, 88);\n"
"}\n"
"QTextEdit:focus {\n"
"    border: 2px solid rgb(91, 101, 124);\n"
"}\n"
"\n"
"QTextEdit:disabled {\n"
"background-color: rgb(93, 83, 91);\n"
"}")
        self.Max_minimize_iter_textEdit.setObjectName("Max_minimize_iter_textEdit")
        self.formLayout_4.setWidget(1, QtWidgets.QFormLayout.FieldRole, self.Max_minimize_iter_textEdit)
        self.formLayoutWidget_5 = QtWidgets.QWidget(self.groupBox_10)
        self.formLayoutWidget_5.setGeometry(QtCore.QRect(290, 30, 261, 81))
        self.formLayoutWidget_5.setObjectName("formLayoutWidget_5")
        self.formLayout_5 = QtWidgets.QFormLayout(self.formLayoutWidget_5)
        self.formLayout_5.setContentsMargins(3, 8, 3, 4)
        self.formLayout_5.setVerticalSpacing(9)
        self.formLayout_5.setObjectName("formLayout_5")
        self.label_43 = QtWidgets.QLabel(self.formLayoutWidget_5)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_43.sizePolicy().hasHeightForWidth())
        self.label_43.setSizePolicy(sizePolicy)
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(10)
        self.label_43.setFont(font)
        self.label_43.setObjectName("label_43")
        self.formLayout_5.setWidget(0, QtWidgets.QFormLayout.LabelRole, self.label_43)
        self.label_44 = QtWidgets.QLabel(self.formLayoutWidget_5)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_44.sizePolicy().hasHeightForWidth())
        self.label_44.setSizePolicy(sizePolicy)
        self.label_44.setMinimumSize(QtCore.QSize(107, 0))
        self.label_44.setMaximumSize(QtCore.QSize(107, 16777215))
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(10)
        self.label_44.setFont(font)
        self.label_44.setObjectName("label_44")
        self.formLayout_5.setWidget(1, QtWidgets.QFormLayout.LabelRole, self.label_44)
        self.equilubrate_checkBox = QtWidgets.QCheckBox(self.formLayoutWidget_5)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.equilubrate_checkBox.sizePolicy().hasHeightForWidth())
        self.equilubrate_checkBox.setSizePolicy(sizePolicy)
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(10)
        self.equilubrate_checkBox.setFont(font)
        self.equilubrate_checkBox.setStyleSheet("QCheckBox {\n"
"    spacing: 5px;\n"
"    margin-left: 15px;\n"
"}\n"
"\n"
"QCheckBox::indicator {\n"
"    width: 15px;\n"
"    height: 15px;\n"
"}\n"
"\n"
"QCheckBox::indicator:unchecked {\n"
"border-image: url(:/16x16/icons/16x16/cil-circle.png);\n"
"\n"
"}\n"
"\n"
"\n"
"QCheckBox::indicator:checked {\n"
"    border-image: url(:/20x20/icons/20x20/cil-check.png);\n"
"}\n"
"\n"
"\n"
"")
        self.equilubrate_checkBox.setText("")
        self.equilubrate_checkBox.setChecked(True)
        self.equilubrate_checkBox.setObjectName("equilubrate_checkBox")
        self.formLayout_5.setWidget(0, QtWidgets.QFormLayout.FieldRole, self.equilubrate_checkBox)
        self.Max_equilubrate_steps_textEdit = QtWidgets.QTextEdit(self.formLayoutWidget_5)
        self.Max_equilubrate_steps_textEdit.setMinimumSize(QtCore.QSize(0, 38))
        self.Max_equilubrate_steps_textEdit.setMaximumSize(QtCore.QSize(16777215, 38))
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(10)
        self.Max_equilubrate_steps_textEdit.setFont(font)
        self.Max_equilubrate_steps_textEdit.setStyleSheet("QTextEdit {\n"
"    background-color: rgb(27, 29, 35);\n"
"    border-radius: 5px;\n"
"    padding: 10px;\n"
"}\n"
"QTextEdit:hover {\n"
"    border: 2px solid rgb(64, 71, 88);\n"
"}\n"
"QTextEdit:focus {\n"
"    border: 2px solid rgb(91, 101, 124);\n"
"}\n"
"\n"
"QTextEdit:disabled {\n"
"background-color: rgb(93, 83, 91);\n"
"}")
        self.Max_equilubrate_steps_textEdit.setObjectName("Max_equilubrate_steps_textEdit")
        self.formLayout_5.setWidget(1, QtWidgets.QFormLayout.FieldRole, self.Max_equilubrate_steps_textEdit)
        self.horizontalLayoutWidget_7 = QtWidgets.QWidget(self.groupBox_8)
        self.horizontalLayoutWidget_7.setGeometry(QtCore.QRect(19, 30, 331, 50))
        self.horizontalLayoutWidget_7.setObjectName("horizontalLayoutWidget_7")
        self.horizontalLayout_14 = QtWidgets.QHBoxLayout(self.horizontalLayoutWidget_7)
        self.horizontalLayout_14.setContentsMargins(3, 8, 3, 4)
        self.horizontalLayout_14.setObjectName("horizontalLayout_14")
        self.label_19 = QtWidgets.QLabel(self.horizontalLayoutWidget_7)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_19.sizePolicy().hasHeightForWidth())
        self.label_19.setSizePolicy(sizePolicy)
        self.label_19.setMaximumSize(QtCore.QSize(107, 16777215))
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(10)
        self.label_19.setFont(font)
        self.label_19.setObjectName("label_19")
        self.horizontalLayout_14.addWidget(self.label_19)
        self.Number_of_steps_spinBox = QtWidgets.QSpinBox(self.horizontalLayoutWidget_7)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.Number_of_steps_spinBox.sizePolicy().hasHeightForWidth())
        self.Number_of_steps_spinBox.setSizePolicy(sizePolicy)
        self.Number_of_steps_spinBox.setMinimumSize(QtCore.QSize(0, 38))
        self.Number_of_steps_spinBox.setMaximumSize(QtCore.QSize(16777215, 38))
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(10)
        self.Number_of_steps_spinBox.setFont(font)
        self.Number_of_steps_spinBox.setStyleSheet("QSpinBox{\n"
"    background-color: rgb(27, 29, 35);\n"
"    border-radius: 5px;\n"
"    border: 2px solid rgb(27, 29, 35);\n"
"    padding: 5px;\n"
"    padding-left: 5px;\n"
"}\n"
"QSpinBox:hover{\n"
"    border: 2px solid rgb(64, 71, 88);\n"
"}\n"
"\n"
"QSpinBox QAbstractItemView {\n"
"    color: rgb(85, 170, 255);    \n"
"    background-color: rgb(27, 29, 35);\n"
"    padding: 5px;\n"
"    selection-background-color: rgb(39, 44, 54);\n"
"}")
        self.Number_of_steps_spinBox.setMaximum(500000000)
        self.Number_of_steps_spinBox.setProperty("value", 1000)
        self.Number_of_steps_spinBox.setObjectName("Number_of_steps_spinBox")
        self.horizontalLayout_14.addWidget(self.Number_of_steps_spinBox)
        self.verticalLayout_11.addWidget(self.groupBox_8)
        self.groupBox_9 = QtWidgets.QGroupBox(self.layoutWidget1)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.groupBox_9.sizePolicy().hasHeightForWidth())
        self.groupBox_9.setSizePolicy(sizePolicy)
        self.groupBox_9.setMinimumSize(QtCore.QSize(0, 500))
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(10)
        self.groupBox_9.setFont(font)
        self.groupBox_9.setStyleSheet("QGroupBox {\n"
"  \n"
"  border: 1px solid rgb(220, 220, 220);\n"
"  border-radius: 10px;\n"
"  padding: 4px;\n"
"  margin-top: 16px;\n"
"}\n"
"\n"
"QGroupBox::title {\n"
"  subcontrol-origin: margin;\n"
"  subcontrol-position: top left;\n"
"  left: 10px;\n"
"  padding-left: 3px;\n"
"  padding-right: 5px;\n"
"  padding-top: 8px;\n"
"  padding-bottom: 16px;\n"
"}\n"
"\n"
"QGroupBox::indicator {\n"
"  margin-left: 2px;\n"
"  height: 16px;\n"
"  width: 16px;\n"
"}")
        self.groupBox_9.setObjectName("groupBox_9")
        self.State_Data_Reporter_Options_groupBox = QtWidgets.QGroupBox(self.groupBox_9)
        self.State_Data_Reporter_Options_groupBox.setGeometry(QtCore.QRect(9, 400, 351, 91))
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(10)
        self.State_Data_Reporter_Options_groupBox.setFont(font)
        self.State_Data_Reporter_Options_groupBox.setStyleSheet("QWidget:disabled {\n"
"  background-color: #19232D;\n"
"  color: #787878;\n"
"  selection-background-color: #14506E;\n"
"  selection-color: #787878;\n"
"}")
        self.State_Data_Reporter_Options_groupBox.setObjectName("State_Data_Reporter_Options_groupBox")
        self.horizontalLayoutWidget_8 = QtWidgets.QWidget(self.State_Data_Reporter_Options_groupBox)
        self.horizontalLayoutWidget_8.setGeometry(QtCore.QRect(10, 30, 321, 46))
        self.horizontalLayoutWidget_8.setObjectName("horizontalLayoutWidget_8")
        self.horizontalLayout_15 = QtWidgets.QHBoxLayout(self.horizontalLayoutWidget_8)
        self.horizontalLayout_15.setContentsMargins(3, 4, 3, 4)
        self.horizontalLayout_15.setObjectName("horizontalLayout_15")
        self.label_26 = QtWidgets.QLabel(self.horizontalLayoutWidget_8)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_26.sizePolicy().hasHeightForWidth())
        self.label_26.setSizePolicy(sizePolicy)
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(10)
        self.label_26.setFont(font)
        self.label_26.setObjectName("label_26")
        self.horizontalLayout_15.addWidget(self.label_26)
        self.StateData_frequency_textEdit = QtWidgets.QTextEdit(self.horizontalLayoutWidget_8)
        self.StateData_frequency_textEdit.setMinimumSize(QtCore.QSize(0, 38))
        self.StateData_frequency_textEdit.setSizeIncrement(QtCore.QSize(0, 38))
        self.StateData_frequency_textEdit.setStyleSheet("QTextEdit {\n"
"    background-color: rgb(27, 29, 35);\n"
"    border-radius: 5px;\n"
"    padding: 10px;\n"
"}\n"
"QTextEdit:hover {\n"
"    border: 2px solid rgb(64, 71, 88);\n"
"}\n"
"QTextEdit:focus {\n"
"    border: 2px solid rgb(91, 101, 124);\n"
"}")
        self.StateData_frequency_textEdit.setObjectName("StateData_frequency_textEdit")
        self.horizontalLayout_15.addWidget(self.StateData_frequency_textEdit)
        self.DCD_Reporter_Options_groupBox = QtWidgets.QGroupBox(self.groupBox_9)
        self.DCD_Reporter_Options_groupBox.setGeometry(QtCore.QRect(10, 110, 511, 131))
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(10)
        self.DCD_Reporter_Options_groupBox.setFont(font)
        self.DCD_Reporter_Options_groupBox.setStyleSheet("QWidget:disabled {\n"
"  background-color: #19232D;\n"
"  color: #787878;\n"
"  selection-background-color: #14506E;\n"
"  selection-color: #787878;\n"
"}")
        self.DCD_Reporter_Options_groupBox.setObjectName("DCD_Reporter_Options_groupBox")
        self.formLayoutWidget_6 = QtWidgets.QWidget(self.DCD_Reporter_Options_groupBox)
        self.formLayoutWidget_6.setGeometry(QtCore.QRect(10, 29, 491, 92))
        self.formLayoutWidget_6.setObjectName("formLayoutWidget_6")
        self.formLayout_6 = QtWidgets.QFormLayout(self.formLayoutWidget_6)
        self.formLayout_6.setContentsMargins(3, 3, 3, 3)
        self.formLayout_6.setObjectName("formLayout_6")
        self.label_24 = QtWidgets.QLabel(self.formLayoutWidget_6)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_24.sizePolicy().hasHeightForWidth())
        self.label_24.setSizePolicy(sizePolicy)
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(10)
        self.label_24.setFont(font)
        self.label_24.setObjectName("label_24")
        self.formLayout_6.setWidget(0, QtWidgets.QFormLayout.LabelRole, self.label_24)
        self.label_25 = QtWidgets.QLabel(self.formLayoutWidget_6)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_25.sizePolicy().hasHeightForWidth())
        self.label_25.setSizePolicy(sizePolicy)
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(10)
        self.label_25.setFont(font)
        self.label_25.setObjectName("label_25")
        self.formLayout_6.setWidget(1, QtWidgets.QFormLayout.LabelRole, self.label_25)
        self.DCD_write_freq_textEdit = QtWidgets.QTextEdit(self.formLayoutWidget_6)
        self.DCD_write_freq_textEdit.setMinimumSize(QtCore.QSize(0, 38))
        self.DCD_write_freq_textEdit.setMaximumSize(QtCore.QSize(16777215, 38))
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(10)
        self.DCD_write_freq_textEdit.setFont(font)
        self.DCD_write_freq_textEdit.setStyleSheet("QTextEdit {\n"
"    background-color: rgb(27, 29, 35);\n"
"    border-radius: 5px;\n"
"    padding: 10px;\n"
"}\n"
"QTextEdit:hover {\n"
"    border: 2px solid rgb(64, 71, 88);\n"
"}\n"
"QTextEdit:focus {\n"
"    border: 2px solid rgb(91, 101, 124);\n"
"}")
        self.DCD_write_freq_textEdit.setObjectName("DCD_write_freq_textEdit")
        self.formLayout_6.setWidget(0, QtWidgets.QFormLayout.FieldRole, self.DCD_write_freq_textEdit)
        self.DCD_Output_Name_textEdit = QtWidgets.QTextEdit(self.formLayoutWidget_6)
        self.DCD_Output_Name_textEdit.setMinimumSize(QtCore.QSize(0, 38))
        self.DCD_Output_Name_textEdit.setMaximumSize(QtCore.QSize(16777215, 38))
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(10)
        self.DCD_Output_Name_textEdit.setFont(font)
        self.DCD_Output_Name_textEdit.setStyleSheet("QTextEdit {\n"
"    background-color: rgb(27, 29, 35);\n"
"    border-radius: 5px;\n"
"    padding: 10px;\n"
"}\n"
"QTextEdit:hover {\n"
"    border: 2px solid rgb(64, 71, 88);\n"
"}\n"
"QTextEdit:focus {\n"
"    border: 2px solid rgb(91, 101, 124);\n"
"}")
        self.DCD_Output_Name_textEdit.setObjectName("DCD_Output_Name_textEdit")
        self.formLayout_6.setWidget(1, QtWidgets.QFormLayout.FieldRole, self.DCD_Output_Name_textEdit)
        self.XTC_Reporter_Options_groupBox = QtWidgets.QGroupBox(self.groupBox_9)
        self.XTC_Reporter_Options_groupBox.setEnabled(False)
        self.XTC_Reporter_Options_groupBox.setGeometry(QtCore.QRect(10, 250, 511, 131))
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(10)
        self.XTC_Reporter_Options_groupBox.setFont(font)
        self.XTC_Reporter_Options_groupBox.setStyleSheet("QWidget:disabled {\n"
"  background-color: #19232D;\n"
"  color: #787878;\n"
"  selection-background-color: #14506E;\n"
"  selection-color: #787878;\n"
"}")
        self.XTC_Reporter_Options_groupBox.setObjectName("XTC_Reporter_Options_groupBox")
        self.formLayoutWidget_10 = QtWidgets.QWidget(self.XTC_Reporter_Options_groupBox)
        self.formLayoutWidget_10.setGeometry(QtCore.QRect(10, 29, 481, 92))
        self.formLayoutWidget_10.setObjectName("formLayoutWidget_10")
        self.formLayout_9 = QtWidgets.QFormLayout(self.formLayoutWidget_10)
        self.formLayout_9.setContentsMargins(3, 3, 3, 3)
        self.formLayout_9.setObjectName("formLayout_9")
        self.label_36 = QtWidgets.QLabel(self.formLayoutWidget_10)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_36.sizePolicy().hasHeightForWidth())
        self.label_36.setSizePolicy(sizePolicy)
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(10)
        self.label_36.setFont(font)
        self.label_36.setObjectName("label_36")
        self.formLayout_9.setWidget(0, QtWidgets.QFormLayout.LabelRole, self.label_36)
        self.label_37 = QtWidgets.QLabel(self.formLayoutWidget_10)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_37.sizePolicy().hasHeightForWidth())
        self.label_37.setSizePolicy(sizePolicy)
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(10)
        self.label_37.setFont(font)
        self.label_37.setObjectName("label_37")
        self.formLayout_9.setWidget(1, QtWidgets.QFormLayout.LabelRole, self.label_37)
        self.XTC_write_freq_textEdit = QtWidgets.QTextEdit(self.formLayoutWidget_10)
        self.XTC_write_freq_textEdit.setMinimumSize(QtCore.QSize(0, 38))
        self.XTC_write_freq_textEdit.setMaximumSize(QtCore.QSize(16777215, 38))
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(10)
        self.XTC_write_freq_textEdit.setFont(font)
        self.XTC_write_freq_textEdit.setStyleSheet("QTextEdit {\n"
"    background-color: rgb(27, 29, 35);\n"
"    border-radius: 5px;\n"
"    padding: 10px;\n"
"}\n"
"QTextEdit:hover {\n"
"    border: 2px solid rgb(64, 71, 88);\n"
"}\n"
"QTextEdit:focus {\n"
"    border: 2px solid rgb(91, 101, 124);\n"
"}")
        self.XTC_write_freq_textEdit.setObjectName("XTC_write_freq_textEdit")
        self.formLayout_9.setWidget(0, QtWidgets.QFormLayout.FieldRole, self.XTC_write_freq_textEdit)
        self.XTC_Output_Name_textEdit = QtWidgets.QTextEdit(self.formLayoutWidget_10)
        self.XTC_Output_Name_textEdit.setMinimumSize(QtCore.QSize(0, 38))
        self.XTC_Output_Name_textEdit.setMaximumSize(QtCore.QSize(16777215, 38))
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(10)
        self.XTC_Output_Name_textEdit.setFont(font)
        self.XTC_Output_Name_textEdit.setStyleSheet("QTextEdit {\n"
"    background-color: rgb(27, 29, 35);\n"
"    border-radius: 5px;\n"
"    padding: 10px;\n"
"}\n"
"QTextEdit:hover {\n"
"    border: 2px solid rgb(64, 71, 88);\n"
"}\n"
"QTextEdit:focus {\n"
"    border: 2px solid rgb(91, 101, 124);\n"
"}")
        self.XTC_Output_Name_textEdit.setObjectName("XTC_Output_Name_textEdit")
        self.formLayout_9.setWidget(1, QtWidgets.QFormLayout.FieldRole, self.XTC_Output_Name_textEdit)
        self.gridLayoutWidget_5 = QtWidgets.QWidget(self.groupBox_9)
        self.gridLayoutWidget_5.setGeometry(QtCore.QRect(20, 30, 351, 71))
        self.gridLayoutWidget_5.setObjectName("gridLayoutWidget_5")
        self.gridLayout_6 = QtWidgets.QGridLayout(self.gridLayoutWidget_5)
        self.gridLayout_6.setContentsMargins(3, 3, 3, 3)
        self.gridLayout_6.setObjectName("gridLayout_6")
        self.DCD_Reporter_checkBox = QtWidgets.QCheckBox(self.gridLayoutWidget_5)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.DCD_Reporter_checkBox.sizePolicy().hasHeightForWidth())
        self.DCD_Reporter_checkBox.setSizePolicy(sizePolicy)
        self.DCD_Reporter_checkBox.setMinimumSize(QtCore.QSize(40, 0))
        self.DCD_Reporter_checkBox.setMaximumSize(QtCore.QSize(40, 16777215))
        self.DCD_Reporter_checkBox.setStyleSheet("QCheckBox {\n"
"    spacing: 5px;\n"
"    margin-left: 15px;\n"
"}\n"
"\n"
"QCheckBox::indicator {\n"
"    width: 15px;\n"
"    height: 15px;\n"
"}\n"
"\n"
"QCheckBox::indicator:unchecked {\n"
"border-image: url(:/16x16/icons/16x16/cil-circle.png);\n"
"\n"
"}\n"
"\n"
"\n"
"QCheckBox::indicator:checked {\n"
"    border-image: url(:/20x20/icons/20x20/cil-check.png);\n"
"}\n"
"\n"
"\n"
"")
        self.DCD_Reporter_checkBox.setText("")
        self.DCD_Reporter_checkBox.setChecked(True)
        self.DCD_Reporter_checkBox.setObjectName("DCD_Reporter_checkBox")
        self.gridLayout_6.addWidget(self.DCD_Reporter_checkBox, 0, 1, 1, 1)
        self.label_22 = QtWidgets.QLabel(self.gridLayoutWidget_5)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_22.sizePolicy().hasHeightForWidth())
        self.label_22.setSizePolicy(sizePolicy)
        self.label_22.setMaximumSize(QtCore.QSize(115, 16777215))
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(10)
        self.label_22.setFont(font)
        self.label_22.setObjectName("label_22")
        self.gridLayout_6.addWidget(self.label_22, 0, 0, 1, 1)
        self.label_40 = QtWidgets.QLabel(self.gridLayoutWidget_5)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_40.sizePolicy().hasHeightForWidth())
        self.label_40.setSizePolicy(sizePolicy)
        self.label_40.setMinimumSize(QtCore.QSize(0, 0))
        self.label_40.setMaximumSize(QtCore.QSize(115, 16777215))
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(10)
        self.label_40.setFont(font)
        self.label_40.setObjectName("label_40")
        self.gridLayout_6.addWidget(self.label_40, 0, 2, 1, 1)
        self.XTC_Reporter_checkBox = QtWidgets.QCheckBox(self.gridLayoutWidget_5)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.XTC_Reporter_checkBox.sizePolicy().hasHeightForWidth())
        self.XTC_Reporter_checkBox.setSizePolicy(sizePolicy)
        self.XTC_Reporter_checkBox.setMaximumSize(QtCore.QSize(40, 16777215))
        self.XTC_Reporter_checkBox.setStyleSheet("QCheckBox {\n"
"    spacing: 5px;\n"
"    margin-left: 15px;\n"
"}\n"
"\n"
"QCheckBox::indicator {\n"
"    width: 15px;\n"
"    height: 15px;\n"
"}\n"
"\n"
"QCheckBox::indicator:unchecked {\n"
"border-image: url(:/16x16/icons/16x16/cil-circle.png);\n"
"\n"
"}\n"
"\n"
"\n"
"QCheckBox::indicator:checked {\n"
"    border-image: url(:/20x20/icons/20x20/cil-check.png);\n"
"}\n"
"\n"
"\n"
"")
        self.XTC_Reporter_checkBox.setText("")
        self.XTC_Reporter_checkBox.setChecked(False)
        self.XTC_Reporter_checkBox.setObjectName("XTC_Reporter_checkBox")
        self.gridLayout_6.addWidget(self.XTC_Reporter_checkBox, 0, 3, 1, 1)
        self.label_23 = QtWidgets.QLabel(self.gridLayoutWidget_5)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_23.sizePolicy().hasHeightForWidth())
        self.label_23.setSizePolicy(sizePolicy)
        self.label_23.setMinimumSize(QtCore.QSize(115, 0))
        self.label_23.setMaximumSize(QtCore.QSize(115, 16777215))
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(10)
        self.label_23.setFont(font)
        self.label_23.setObjectName("label_23")
        self.gridLayout_6.addWidget(self.label_23, 1, 0, 1, 1)
        self.State_Data_Reporter_checkBox = QtWidgets.QCheckBox(self.gridLayoutWidget_5)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.State_Data_Reporter_checkBox.sizePolicy().hasHeightForWidth())
        self.State_Data_Reporter_checkBox.setSizePolicy(sizePolicy)
        self.State_Data_Reporter_checkBox.setMinimumSize(QtCore.QSize(85, 0))
        self.State_Data_Reporter_checkBox.setMaximumSize(QtCore.QSize(16777215, 16777215))
        self.State_Data_Reporter_checkBox.setStyleSheet("QCheckBox {\n"
"    spacing: 5px;\n"
"    margin-left: 15px;\n"
"}\n"
"\n"
"QCheckBox::indicator {\n"
"    width: 15px;\n"
"    height: 15px;\n"
"}\n"
"\n"
"QCheckBox::indicator:unchecked {\n"
"border-image: url(:/16x16/icons/16x16/cil-circle.png);\n"
"\n"
"}\n"
"\n"
"\n"
"QCheckBox::indicator:checked {\n"
"    border-image: url(:/20x20/icons/20x20/cil-check.png);\n"
"}\n"
"\n"
"\n"
"")
        self.State_Data_Reporter_checkBox.setText("")
        self.State_Data_Reporter_checkBox.setChecked(True)
        self.State_Data_Reporter_checkBox.setObjectName("State_Data_Reporter_checkBox")
        self.gridLayout_6.addWidget(self.State_Data_Reporter_checkBox, 1, 1, 1, 1)
        self.verticalLayout_11.addWidget(self.groupBox_9)
        self.tabWidget.addTab(self.tab_4, "")
        self.verticalLayout_21.addWidget(self.tabWidget)
        self.Perturbation_TabWidget.addTab(self.Perturb_Advanced, "")
        self.horizontalLayout_12.addWidget(self.Perturbation_TabWidget)
        self.show_navigation = QtWidgets.QPushButton(self.page_home)
        self.show_navigation.setEnabled(True)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.show_navigation.sizePolicy().hasHeightForWidth())
        self.show_navigation.setSizePolicy(sizePolicy)
        self.show_navigation.setMaximumSize(QtCore.QSize(20, 61))
        self.show_navigation.setStyleSheet(" QPushButton \n"
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
        self.show_navigation.setText("")
        icon10 = QtGui.QIcon()
        icon10.addPixmap(QtGui.QPixmap(":/16x16/icons/16x16/cil-chevron-circle-right-alt.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.show_navigation.setIcon(icon10)
        self.show_navigation.setObjectName("show_navigation")
        self.horizontalLayout_12.addWidget(self.show_navigation)
        self.Pymol_Widget = QtWidgets.QFrame(self.page_home)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.Pymol_Widget.sizePolicy().hasHeightForWidth())
        self.Pymol_Widget.setSizePolicy(sizePolicy)
        self.Pymol_Widget.setMinimumSize(QtCore.QSize(0, 0))
        self.Pymol_Widget.setStyleSheet("QFrame\n"
"{\n"
"    border: 1px solid black;\n"
"    border-radius: 5px;\n"
"    border-top-color: rgb(157, 90, 198);\n"
"    border-left-color: rgb(157, 90, 198);\n"
"    border-bottom-color: rgb(157, 90, 198);\n"
"    border-right-color: rgb(157, 90, 198);\n"
"    margin-top: 25px\n"
" \n"
"}\n"
"")
        self.Pymol_Widget.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.Pymol_Widget.setFrameShadow(QtWidgets.QFrame.Raised)
        self.Pymol_Widget.setObjectName("Pymol_Widget")
        self.horizontalLayout_12.addWidget(self.Pymol_Widget)
        self.hide_navigation = QtWidgets.QPushButton(self.page_home)
        self.hide_navigation.setMaximumSize(QtCore.QSize(20, 61))
        self.hide_navigation.setStyleSheet(" QPushButton \n"
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
        self.hide_navigation.setText("")
        self.hide_navigation.setIcon(icon10)
        self.hide_navigation.setObjectName("hide_navigation")
        self.horizontalLayout_12.addWidget(self.hide_navigation)
        self.visualization_settings_groupBox = QtWidgets.QGroupBox(self.page_home)
        self.visualization_settings_groupBox.setMinimumSize(QtCore.QSize(170, 0))
        self.visualization_settings_groupBox.setMaximumSize(QtCore.QSize(170, 16777215))
        self.visualization_settings_groupBox.setStyleSheet("QGroupBox\n"
"{\n"
"    border: 1px solid black;\n"
"    border-radius: 5px;\n"
"    border-top-color: rgb(157, 90, 198);\n"
"    border-left-color: rgb(157, 90, 198);\n"
"    border-bottom-color: rgb(157, 90, 198);\n"
"    border-right-color: rgb(157, 90, 198);\n"
"    margin-top: 25px\n"
" \n"
"}")
        self.visualization_settings_groupBox.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignVCenter)
        self.visualization_settings_groupBox.setObjectName("visualization_settings_groupBox")
        self.gridLayoutWidget_6 = QtWidgets.QWidget(self.visualization_settings_groupBox)
        self.gridLayoutWidget_6.setGeometry(QtCore.QRect(11, 50, 151, 251))
        self.gridLayoutWidget_6.setObjectName("gridLayoutWidget_6")
        self.verticalLayout_19 = QtWidgets.QVBoxLayout(self.gridLayoutWidget_6)
        self.verticalLayout_19.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout_19.setObjectName("verticalLayout_19")
        self.activate_pymol_navigation = QtWidgets.QPushButton(self.gridLayoutWidget_6)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.activate_pymol_navigation.sizePolicy().hasHeightForWidth())
        self.activate_pymol_navigation.setSizePolicy(sizePolicy)
        self.activate_pymol_navigation.setMinimumSize(QtCore.QSize(149, 33))
        self.activate_pymol_navigation.setMaximumSize(QtCore.QSize(110, 33))
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(9)
        self.activate_pymol_navigation.setFont(font)
        self.activate_pymol_navigation.setStyleSheet("            QPushButton \n"
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
        icon11 = QtGui.QIcon()
        icon11.addPixmap(QtGui.QPixmap(":/24x24/icons/24x24/cil-cursor.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.activate_pymol_navigation.setIcon(icon11)
        self.activate_pymol_navigation.setObjectName("activate_pymol_navigation")
        self.verticalLayout_19.addWidget(self.activate_pymol_navigation)
        self.deactivate_pymol_navigation = QtWidgets.QPushButton(self.gridLayoutWidget_6)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.deactivate_pymol_navigation.sizePolicy().hasHeightForWidth())
        self.deactivate_pymol_navigation.setSizePolicy(sizePolicy)
        self.deactivate_pymol_navigation.setMinimumSize(QtCore.QSize(149, 33))
        self.deactivate_pymol_navigation.setMaximumSize(QtCore.QSize(110, 33))
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(9)
        self.deactivate_pymol_navigation.setFont(font)
        self.deactivate_pymol_navigation.setStyleSheet("            QPushButton \n"
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
        icon12.addPixmap(QtGui.QPixmap(":/20x20/icons/20x20/cil-x.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.deactivate_pymol_navigation.setIcon(icon12)
        self.deactivate_pymol_navigation.setObjectName("deactivate_pymol_navigation")
        self.verticalLayout_19.addWidget(self.deactivate_pymol_navigation)
        self.refresh_pushButton = QtWidgets.QPushButton(self.gridLayoutWidget_6)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.refresh_pushButton.sizePolicy().hasHeightForWidth())
        self.refresh_pushButton.setSizePolicy(sizePolicy)
        self.refresh_pushButton.setMinimumSize(QtCore.QSize(149, 33))
        self.refresh_pushButton.setMaximumSize(QtCore.QSize(110, 33))
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(9)
        self.refresh_pushButton.setFont(font)
        self.refresh_pushButton.setStyleSheet("            QPushButton \n"
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
        icon13.addPixmap(QtGui.QPixmap(":/16x16/icons/16x16/cil-reload.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.refresh_pushButton.setIcon(icon13)
        self.refresh_pushButton.setObjectName("refresh_pushButton")
        self.verticalLayout_19.addWidget(self.refresh_pushButton)
        self.ss_beatiful_snapshoot = QtWidgets.QPushButton(self.gridLayoutWidget_6)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.ss_beatiful_snapshoot.sizePolicy().hasHeightForWidth())
        self.ss_beatiful_snapshoot.setSizePolicy(sizePolicy)
        self.ss_beatiful_snapshoot.setMinimumSize(QtCore.QSize(149, 33))
        self.ss_beatiful_snapshoot.setMaximumSize(QtCore.QSize(110, 33))
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(9)
        self.ss_beatiful_snapshoot.setFont(font)
        self.ss_beatiful_snapshoot.setStyleSheet("            QPushButton \n"
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
        icon14.addPixmap(QtGui.QPixmap(":/16x16/icons/16x16/cil-camera.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.ss_beatiful_snapshoot.setIcon(icon14)
        self.ss_beatiful_snapshoot.setObjectName("ss_beatiful_snapshoot")
        self.verticalLayout_19.addWidget(self.ss_beatiful_snapshoot)
        self.save_as_png_pushButton = QtWidgets.QPushButton(self.gridLayoutWidget_6)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.save_as_png_pushButton.sizePolicy().hasHeightForWidth())
        self.save_as_png_pushButton.setSizePolicy(sizePolicy)
        self.save_as_png_pushButton.setMinimumSize(QtCore.QSize(149, 33))
        self.save_as_png_pushButton.setMaximumSize(QtCore.QSize(110, 33))
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(9)
        self.save_as_png_pushButton.setFont(font)
        self.save_as_png_pushButton.setStyleSheet("            QPushButton \n"
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
        icon15.addPixmap(QtGui.QPixmap(":/16x16/icons/16x16/cil-data-transfer-down.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.save_as_png_pushButton.setIcon(icon15)
        self.save_as_png_pushButton.setObjectName("save_as_png_pushButton")
        self.verticalLayout_19.addWidget(self.save_as_png_pushButton)
        self.figure_settings_groupBox = QtWidgets.QGroupBox(self.visualization_settings_groupBox)
        self.figure_settings_groupBox.setGeometry(QtCore.QRect(-1, 300, 171, 401))
        self.figure_settings_groupBox.setObjectName("figure_settings_groupBox")
        self.verticalLayoutWidget_2 = QtWidgets.QWidget(self.figure_settings_groupBox)
        self.verticalLayoutWidget_2.setGeometry(QtCore.QRect(10, 54, 151, 311))
        self.verticalLayoutWidget_2.setObjectName("verticalLayoutWidget_2")
        self.figure_settings_verticalLayout = QtWidgets.QVBoxLayout(self.verticalLayoutWidget_2)
        self.figure_settings_verticalLayout.setContentsMargins(0, 0, 0, 0)
        self.figure_settings_verticalLayout.setObjectName("figure_settings_verticalLayout")
        self.pymol_width_label = QtWidgets.QLabel(self.verticalLayoutWidget_2)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.pymol_width_label.sizePolicy().hasHeightForWidth())
        self.pymol_width_label.setSizePolicy(sizePolicy)
        self.pymol_width_label.setMinimumSize(QtCore.QSize(50, 22))
        self.pymol_width_label.setMaximumSize(QtCore.QSize(16777215, 22))
        self.pymol_width_label.setStyleSheet("QLabel {\n"
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
        self.pymol_width_label.setObjectName("pymol_width_label")
        self.figure_settings_verticalLayout.addWidget(self.pymol_width_label)
        self.width_horizontalSlider = QtWidgets.QSlider(self.verticalLayoutWidget_2)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.width_horizontalSlider.sizePolicy().hasHeightForWidth())
        self.width_horizontalSlider.setSizePolicy(sizePolicy)
        self.width_horizontalSlider.setStyleSheet("\n"
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
        self.width_horizontalSlider.setMinimum(800)
        self.width_horizontalSlider.setMaximum(1920)
        self.width_horizontalSlider.setSingleStep(20)
        self.width_horizontalSlider.setPageStep(20)
        self.width_horizontalSlider.setProperty("value", 1200)
        self.width_horizontalSlider.setOrientation(QtCore.Qt.Horizontal)
        self.width_horizontalSlider.setInvertedAppearance(False)
        self.width_horizontalSlider.setInvertedControls(False)
        self.width_horizontalSlider.setTickPosition(QtWidgets.QSlider.NoTicks)
        self.width_horizontalSlider.setTickInterval(20)
        self.width_horizontalSlider.setObjectName("width_horizontalSlider")
        self.figure_settings_verticalLayout.addWidget(self.width_horizontalSlider)
        self.pymol_height_label = QtWidgets.QLabel(self.verticalLayoutWidget_2)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.pymol_height_label.sizePolicy().hasHeightForWidth())
        self.pymol_height_label.setSizePolicy(sizePolicy)
        self.pymol_height_label.setMinimumSize(QtCore.QSize(50, 22))
        self.pymol_height_label.setMaximumSize(QtCore.QSize(16777215, 22))
        self.pymol_height_label.setStyleSheet("QLabel {\n"
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
        self.pymol_height_label.setObjectName("pymol_height_label")
        self.figure_settings_verticalLayout.addWidget(self.pymol_height_label)
        self.height_horizontalSlider = QtWidgets.QSlider(self.verticalLayoutWidget_2)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.height_horizontalSlider.sizePolicy().hasHeightForWidth())
        self.height_horizontalSlider.setSizePolicy(sizePolicy)
        self.height_horizontalSlider.setStyleSheet("QSlider::handle:horizontal \n"
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
        self.height_horizontalSlider.setMinimum(600)
        self.height_horizontalSlider.setMaximum(1080)
        self.height_horizontalSlider.setSingleStep(20)
        self.height_horizontalSlider.setPageStep(20)
        self.height_horizontalSlider.setProperty("value", 1080)
        self.height_horizontalSlider.setOrientation(QtCore.Qt.Horizontal)
        self.height_horizontalSlider.setTickInterval(20)
        self.height_horizontalSlider.setObjectName("height_horizontalSlider")
        self.figure_settings_verticalLayout.addWidget(self.height_horizontalSlider)
        self.pymol_dpi_label = QtWidgets.QLabel(self.verticalLayoutWidget_2)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.pymol_dpi_label.sizePolicy().hasHeightForWidth())
        self.pymol_dpi_label.setSizePolicy(sizePolicy)
        self.pymol_dpi_label.setMinimumSize(QtCore.QSize(50, 22))
        self.pymol_dpi_label.setMaximumSize(QtCore.QSize(16777215, 22))
        self.pymol_dpi_label.setStyleSheet("QLabel {\n"
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
        self.pymol_dpi_label.setObjectName("pymol_dpi_label")
        self.figure_settings_verticalLayout.addWidget(self.pymol_dpi_label)
        self.dpi_horizontalSlider = QtWidgets.QSlider(self.verticalLayoutWidget_2)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.dpi_horizontalSlider.sizePolicy().hasHeightForWidth())
        self.dpi_horizontalSlider.setSizePolicy(sizePolicy)
        self.dpi_horizontalSlider.setStyleSheet("QSlider::handle:horizontal \n"
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
        self.dpi_horizontalSlider.setMinimum(100)
        self.dpi_horizontalSlider.setMaximum(300)
        self.dpi_horizontalSlider.setSingleStep(50)
        self.dpi_horizontalSlider.setPageStep(50)
        self.dpi_horizontalSlider.setProperty("value", 150)
        self.dpi_horizontalSlider.setOrientation(QtCore.Qt.Horizontal)
        self.dpi_horizontalSlider.setTickInterval(50)
        self.dpi_horizontalSlider.setObjectName("dpi_horizontalSlider")
        self.figure_settings_verticalLayout.addWidget(self.dpi_horizontalSlider)
        self.pymol_ray_label = QtWidgets.QLabel(self.verticalLayoutWidget_2)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.pymol_ray_label.sizePolicy().hasHeightForWidth())
        self.pymol_ray_label.setSizePolicy(sizePolicy)
        self.pymol_ray_label.setMinimumSize(QtCore.QSize(0, 22))
        self.pymol_ray_label.setMaximumSize(QtCore.QSize(16777215, 22))
        self.pymol_ray_label.setStyleSheet("QLabel {\n"
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
        self.pymol_ray_label.setObjectName("pymol_ray_label")
        self.figure_settings_verticalLayout.addWidget(self.pymol_ray_label)
        self.ray_horizontalSlider = QtWidgets.QSlider(self.verticalLayoutWidget_2)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.ray_horizontalSlider.sizePolicy().hasHeightForWidth())
        self.ray_horizontalSlider.setSizePolicy(sizePolicy)
        self.ray_horizontalSlider.setStyleSheet("QSlider::handle:horizontal \n"
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
        self.ray_horizontalSlider.setMinimum(1)
        self.ray_horizontalSlider.setMaximum(5)
        self.ray_horizontalSlider.setPageStep(1)
        self.ray_horizontalSlider.setProperty("value", 1)
        self.ray_horizontalSlider.setOrientation(QtCore.Qt.Horizontal)
        self.ray_horizontalSlider.setTickInterval(1)
        self.ray_horizontalSlider.setObjectName("ray_horizontalSlider")
        self.figure_settings_verticalLayout.addWidget(self.ray_horizontalSlider)
        self.get_figure_pushButton = QtWidgets.QPushButton(self.verticalLayoutWidget_2)
        self.get_figure_pushButton.setMinimumSize(QtCore.QSize(0, 33))
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(9)
        self.get_figure_pushButton.setFont(font)
        self.get_figure_pushButton.setStyleSheet("QPushButton \n"
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
        icon16 = QtGui.QIcon()
        icon16.addPixmap(QtGui.QPixmap(":/16x16/icons/16x16/cil-image-plus.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.get_figure_pushButton.setIcon(icon16)
        self.get_figure_pushButton.setObjectName("get_figure_pushButton")
        self.figure_settings_verticalLayout.addWidget(self.get_figure_pushButton)
        self.hide_figure_settings_pushButton = QtWidgets.QPushButton(self.figure_settings_groupBox)
        self.hide_figure_settings_pushButton.setEnabled(True)
        self.hide_figure_settings_pushButton.setGeometry(QtCore.QRect(55, 370, 65, 20))
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.hide_figure_settings_pushButton.sizePolicy().hasHeightForWidth())
        self.hide_figure_settings_pushButton.setSizePolicy(sizePolicy)
        self.hide_figure_settings_pushButton.setMinimumSize(QtCore.QSize(65, 20))
        self.hide_figure_settings_pushButton.setMaximumSize(QtCore.QSize(65, 20))
        font = QtGui.QFont()
        font.setWeight(50)
        font.setBold(False)
        self.hide_figure_settings_pushButton.setFont(font)
        self.hide_figure_settings_pushButton.setStyleSheet(" QPushButton \n"
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
        self.hide_figure_settings_pushButton.setText("")
        icon17 = QtGui.QIcon()
        icon17.addPixmap(QtGui.QPixmap(":/16x16/icons/16x16/cil-chevron-circle-up-alt.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.hide_figure_settings_pushButton.setIcon(icon17)
        self.hide_figure_settings_pushButton.setFlat(False)
        self.hide_figure_settings_pushButton.setObjectName("hide_figure_settings_pushButton")
        self.horizontalLayout_12.addWidget(self.visualization_settings_groupBox)
        self.verticalLayout_10.addLayout(self.horizontalLayout_12)
        self.stackedWidget.addWidget(self.page_home)
        self.Simulation_Graphs = QtWidgets.QWidget()
        self.Simulation_Graphs.setObjectName("Simulation_Graphs")
        self.gridLayout_9 = QtWidgets.QGridLayout(self.Simulation_Graphs)
        self.gridLayout_9.setObjectName("gridLayout_9")
        self.graphics_GroupBox = QtWidgets.QGroupBox(self.Simulation_Graphs)
        self.graphics_GroupBox.setObjectName("graphics_GroupBox")
        self.verticalLayout_16 = QtWidgets.QVBoxLayout(self.graphics_GroupBox)
        self.verticalLayout_16.setObjectName("verticalLayout_16")
        self.verticalLayout_22 = QtWidgets.QVBoxLayout()
        self.verticalLayout_22.setObjectName("verticalLayout_22")
        self.verticalLayout_16.addLayout(self.verticalLayout_22)
        self.groupBox_2 = QtWidgets.QGroupBox(self.graphics_GroupBox)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.groupBox_2.sizePolicy().hasHeightForWidth())
        self.groupBox_2.setSizePolicy(sizePolicy)
        self.groupBox_2.setMinimumSize(QtCore.QSize(0, 70))
        self.groupBox_2.setMaximumSize(QtCore.QSize(16777215, 70))
        font = QtGui.QFont()
        font.setFamily("Segoe UI Semibold")
        font.setPointSize(10)
        self.groupBox_2.setFont(font)
        self.groupBox_2.setStyleSheet("QGroupBox {\n"
"  \n"
"  border: 1px solid rgb(220, 220, 220);\n"
"  border-radius: 10px;\n"
"  padding: 2px;\n"
"  margin-top: 16px;\n"
"}\n"
"\n"
"QGroupBox::title {\n"
"  subcontrol-origin: margin;\n"
"  subcontrol-position: top left;\n"
"  left: 10px;\n"
"  padding-left: 3px;\n"
"  padding-right: 5px;\n"
"  padding-top: 8px;\n"
"  padding-bottom: 16px;\n"
"}\n"
"\n"
"QGroupBox::indicator {\n"
"  margin-left: 2px;\n"
"  height: 16px;\n"
"  width: 16px;\n"
"}")
        self.groupBox_2.setObjectName("groupBox_2")
        self.verticalLayout_15 = QtWidgets.QVBoxLayout(self.groupBox_2)
        self.verticalLayout_15.setObjectName("verticalLayout_15")
        self.formLayout_7 = QtWidgets.QFormLayout()
        self.formLayout_7.setContentsMargins(1, 3, 1, 3)
        self.formLayout_7.setObjectName("formLayout_7")
        self.progressBar_decomp = QtWidgets.QProgressBar(self.groupBox_2)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.progressBar_decomp.sizePolicy().hasHeightForWidth())
        self.progressBar_decomp.setSizePolicy(sizePolicy)
        self.progressBar_decomp.setMaximumSize(QtCore.QSize(16777215, 35))
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(10)
        self.progressBar_decomp.setFont(font)
        self.progressBar_decomp.setStyleSheet("QProgressBar {\n"
"    border: 2px solid grey;\n"
"    border-radius: 5px;\n"
"    text-align: center;\n"
"}\n"
"\n"
"QProgressBar::chunk \n"
"{\n"
"      border-radius: 5px;\n"
"    background: #007AA5;\n"
"    margin: 3px;\n"
"\n"
"\n"
"}")
        self.progressBar_decomp.setProperty("value", 0)
        self.progressBar_decomp.setObjectName("progressBar_decomp")
        self.formLayout_7.setWidget(0, QtWidgets.QFormLayout.FieldRole, self.progressBar_decomp)
        self.verticalLayout_15.addLayout(self.formLayout_7)
        self.verticalLayout_16.addWidget(self.groupBox_2)
        self.gridLayout_9.addWidget(self.graphics_GroupBox, 0, 0, 1, 1)
        self.stackedWidget.addWidget(self.Simulation_Graphs)
        self.Analysis = QtWidgets.QWidget()
        self.Analysis.setObjectName("Analysis")
        self.horizontalLayout_33 = QtWidgets.QHBoxLayout(self.Analysis)
        self.horizontalLayout_33.setObjectName("horizontalLayout_33")
        self.analysis_TabWidget = QtWidgets.QTabWidget(self.Analysis)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.analysis_TabWidget.sizePolicy().hasHeightForWidth())
        self.analysis_TabWidget.setSizePolicy(sizePolicy)
        self.analysis_TabWidget.setMinimumSize(QtCore.QSize(600, 0))
        self.analysis_TabWidget.setMaximumSize(QtCore.QSize(16777215, 16777215))
        font = QtGui.QFont()
        font.setFamily("Segoe UI Semibold")
        font.setPointSize(8)
        self.analysis_TabWidget.setFont(font)
        self.analysis_TabWidget.setStyleSheet("\n"
"\n"
"QTabWidget::pane {\n"
"    border: 1px solid black;\n"
"    background: white;\n"
"    background-color: rgb(44, 49, 60);\n"
"    border-radius: 5px;\n"
"    border-top-color: rgb(157, 90, 198);\n"
"    border-right-color: rgb(157, 90, 198);\n"
"    border-bottom-color: rgb(157, 90, 198);\n"
"    border-left-color: rgb(157, 90, 198);\n"
"}\n"
"\n"
"QTabWidget::tab-bar:top {\n"
"    top: 0px;\n"
"\n"
"}\n"
"\n"
"QTabWidget::tab-bar:bottom {\n"
"    bottom: 0px;\n"
"}\n"
"\n"
"QTabWidget::tab-bar:left {\n"
"    right: 0px;\n"
"}\n"
"\n"
"QTabWidget::tab-bar:right {\n"
"    left: 0px;\n"
"}\n"
"\n"
"QTabBar::tab {\n"
"    border: 2px;\n"
"    \n"
"    border-color: rgb(0, 0, 0);\n"
"    border-radius: 2px;\n"
"    background-color: rgb(255, 17, 100);\n"
"\n"
"    \n"
"}\n"
"\n"
"\n"
"QTabBar::tab:selected {\n"
"    background: white;\n"
"    background-color: rgb(200, 5, 100);\n"
"    \n"
"    \n"
"}\n"
"\n"
"QTabBar::tab:!selected {\n"
"    background: white;\n"
"    \n"
"    background-color: rgb(141, 2, 72);\n"
"\n"
"    \n"
"}\n"
"\n"
"QTabBar::tab:!selected:hover {\n"
"    background: #999;\n"
"}\n"
"\n"
"QTabBar::tab:top:!selected {\n"
"    margin-top: 3px;\n"
"}\n"
"\n"
"QTabBar::tab:bottom:!selected {\n"
"    margin-bottom: 3px;\n"
"}\n"
"\n"
"QTabBar::tab:top, QTabBar::tab:bottom {\n"
"    min-width: 8ex;\n"
"    margin-right: -1px;\n"
"    padding: 5px 10px 5px 10px;\n"
"}\n"
"\n"
"QTabBar::tab:top:selected {\n"
"    border-bottom-color: none;\n"
"}\n"
"\n"
"QTabBar::tab:bottom:selected {\n"
"    border-top-color: none;\n"
"}\n"
"\n"
"QTabBar::tab:top:last, QTabBar::tab:bottom:last,\n"
"QTabBar::tab:top:only-one, QTabBar::tab:bottom:only-one {\n"
"    margin-right: 0;\n"
"}\n"
"\n"
"QTabBar::tab:left:!selected {\n"
"    margin-right: 3px;\n"
"}\n"
"\n"
"QTabBar::tab:right:!selected {\n"
"    margin-left: 3px;\n"
"}\n"
"\n"
"QTabBar::tab:left, QTabBar::tab:right {\n"
"    min-height: 8ex;\n"
"    margin-bottom: 1px;\n"
"    padding: 10px 5px 10px 5px;\n"
"}\n"
"\n"
"QTabBar::tab:left:selected {\n"
"    border-left-color: black;\n"
"}\n"
"\n"
"QTabBar::tab:right:selected {\n"
"    border-right-color: black;\n"
"}\n"
"\n"
"QTabBar::tab:left:last, QTabBar::tab:right:last,\n"
"QTabBar::tab:left:only-one, QTabBar::tab:right:only-one {\n"
"    margin-bottom: 0;\n"
"}")
        self.analysis_TabWidget.setObjectName("analysis_TabWidget")
        self.Perturb_Quick_3 = QtWidgets.QWidget()
        self.Perturb_Quick_3.setObjectName("Perturb_Quick_3")
        self.horizontalLayout_34 = QtWidgets.QHBoxLayout(self.Perturb_Quick_3)
        self.horizontalLayout_34.setObjectName("horizontalLayout_34")
        self.formLayout_2 = QtWidgets.QFormLayout()
        self.formLayout_2.setSizeConstraint(QtWidgets.QLayout.SetDefaultConstraint)
        self.formLayout_2.setContentsMargins(-1, 0, -1, -1)
        self.formLayout_2.setVerticalSpacing(0)
        self.formLayout_2.setObjectName("formLayout_2")
        self.groupBox_29 = QtWidgets.QGroupBox(self.Perturb_Quick_3)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.groupBox_29.sizePolicy().hasHeightForWidth())
        self.groupBox_29.setSizePolicy(sizePolicy)
        self.groupBox_29.setMinimumSize(QtCore.QSize(610, 0))
        self.groupBox_29.setMaximumSize(QtCore.QSize(690, 90))
        font = QtGui.QFont()
        font.setFamily("Segoe UI Semibold")
        font.setPointSize(10)
        self.groupBox_29.setFont(font)
        self.groupBox_29.setStyleSheet("QGroupBox {\n"
"  \n"
"  border: 1px solid rgb(220, 220, 220);\n"
"  border-radius: 10px;\n"
"  padding: 4px;\n"
"  margin-top: 16px;\n"
"}\n"
"\n"
"QGroupBox::title {\n"
"  subcontrol-origin: margin;\n"
"  subcontrol-position: top left;\n"
"  left: 10px;\n"
"  padding-left: 3px;\n"
"  padding-right: 5px;\n"
"  padding-top: 8px;\n"
"  padding-bottom: 16px;\n"
"}\n"
"\n"
"QGroupBox::indicator {\n"
"  margin-left: 2px;\n"
"  height: 16px;\n"
"  width: 16px;\n"
"}")
        self.groupBox_29.setObjectName("groupBox_29")
        self.horizontalLayout_38 = QtWidgets.QHBoxLayout(self.groupBox_29)
        self.horizontalLayout_38.setObjectName("horizontalLayout_38")
        self.horizontalLayout_30 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_30.setObjectName("horizontalLayout_30")
        self.stop_pushButton_3 = QtWidgets.QPushButton(self.groupBox_29)
        self.stop_pushButton_3.setMinimumSize(QtCore.QSize(170, 33))
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(9)
        self.stop_pushButton_3.setFont(font)
        self.stop_pushButton_3.setStyleSheet("QPushButton {\n"
"    border: 2px solid rgb(52, 59, 72);\n"
"    border-radius: 5px;    \n"
"    background-color: rgb(0, 180, 95);\n"
"}\n"
"QPushButton:hover {\n"
"    background-color: rgb(0, 152, 84);\n"
"    border: 2px solid rgb(61, 70, 86);\n"
"}\n"
"QPushButton:pressed {    \n"
"    background-color: rgb(220, 5, 100);\n"
"    border: 2px solid rgb(43, 50, 61);\n"
"}")
        self.stop_pushButton_3.setIcon(icon7)
        self.stop_pushButton_3.setObjectName("stop_pushButton_3")
        self.horizontalLayout_30.addWidget(self.stop_pushButton_3)
        self.network_calculate_pushButton = QtWidgets.QPushButton(self.groupBox_29)
        self.network_calculate_pushButton.setMinimumSize(QtCore.QSize(170, 33))
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(9)
        self.network_calculate_pushButton.setFont(font)
        self.network_calculate_pushButton.setStyleSheet("QPushButton {\n"
"    border: 2px solid rgb(52, 59, 72);\n"
"    border-radius: 5px;    \n"
"    background-color: rgb(0, 180, 95);\n"
"}\n"
"QPushButton:hover {\n"
"    background-color: rgb(220, 5, 100);\n"
"    border: 2px solid rgb(61, 70, 86);\n"
"}\n"
"QPushButton:pressed {    \n"
"    background-color: rgb(190, 5, 100);\n"
"    border: 2px solid rgb(43, 50, 61);\n"
"}\n"
"\n"
"")
        icon18 = QtGui.QIcon()
        icon18.addPixmap(QtGui.QPixmap(":/16x16/icons/16x16/cil-arrow-right.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.network_calculate_pushButton.setIcon(icon18)
        self.network_calculate_pushButton.setObjectName("network_calculate_pushButton")
        self.horizontalLayout_30.addWidget(self.network_calculate_pushButton)
        self.horizontalLayout_38.addLayout(self.horizontalLayout_30)
        self.horizontalLayout_28 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_28.setContentsMargins(3, -1, -1, -1)
        self.horizontalLayout_28.setObjectName("horizontalLayout_28")
        self.label_100 = QtWidgets.QLabel(self.groupBox_29)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_100.sizePolicy().hasHeightForWidth())
        self.label_100.setSizePolicy(sizePolicy)
        self.label_100.setMaximumSize(QtCore.QSize(70, 16777215))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.label_100.setFont(font)
        self.label_100.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignVCenter)
        self.label_100.setObjectName("label_100")
        self.horizontalLayout_28.addWidget(self.label_100)
        self.Number_of_thread_for_network_spinBox = QtWidgets.QSpinBox(self.groupBox_29)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.Number_of_thread_for_network_spinBox.sizePolicy().hasHeightForWidth())
        self.Number_of_thread_for_network_spinBox.setSizePolicy(sizePolicy)
        self.Number_of_thread_for_network_spinBox.setMinimumSize(QtCore.QSize(60, 33))
        self.Number_of_thread_for_network_spinBox.setMaximumSize(QtCore.QSize(16777215, 33))
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(10)
        self.Number_of_thread_for_network_spinBox.setFont(font)
        self.Number_of_thread_for_network_spinBox.setLayoutDirection(QtCore.Qt.RightToLeft)
        self.Number_of_thread_for_network_spinBox.setStyleSheet("QSpinBox{\n"
"    background-color: rgb(27, 29, 35);\n"
"    border-radius: 5px;\n"
"    border: 2px solid rgb(27, 29, 35);\n"
"    padding: 5px;\n"
"    padding-left: 5px;\n"
"}\n"
"QSpinBox:hover{\n"
"    border: 2px solid rgb(64, 71, 88);\n"
"}\n"
"QSpinBox QAbstractItemView {\n"
"    color: rgb(85, 170, 255);    \n"
"    background-color: rgb(27, 29, 35);\n"
"    padding: 5px;\n"
"    selection-background-color: rgb(39, 44, 54);\n"
"}")
        self.Number_of_thread_for_network_spinBox.setMinimum(2)
        self.Number_of_thread_for_network_spinBox.setObjectName("Number_of_thread_for_network_spinBox")
        self.horizontalLayout_28.addWidget(self.Number_of_thread_for_network_spinBox)
        self.checkBox_5 = QtWidgets.QCheckBox(self.groupBox_29)
        self.checkBox_5.setMinimumSize(QtCore.QSize(90, 0))
        self.checkBox_5.setMaximumSize(QtCore.QSize(115, 16777215))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.checkBox_5.setFont(font)
        self.checkBox_5.setStyleSheet("QCheckBox {\n"
"    spacing: 5px;\n"
"    margin-left: 15px;\n"
"}\n"
"\n"
"QCheckBox::indicator {\n"
"    width: 15px;\n"
"    height: 15px;\n"
"}\n"
"\n"
"QCheckBox::indicator:unchecked {\n"
"border-image: url(:/16x16/icons/16x16/cil-circle.png);\n"
"\n"
"}\n"
"\n"
"\n"
"QCheckBox::indicator:checked {\n"
"    border-image: url(:/20x20/icons/20x20/cil-check.png);\n"
"}")
        self.checkBox_5.setObjectName("checkBox_5")
        self.horizontalLayout_28.addWidget(self.checkBox_5)
        self.horizontalLayout_38.addLayout(self.horizontalLayout_28)
        self.formLayout_2.setWidget(4, QtWidgets.QFormLayout.SpanningRole, self.groupBox_29)
        self.groupBox_31 = QtWidgets.QGroupBox(self.Perturb_Quick_3)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.groupBox_31.sizePolicy().hasHeightForWidth())
        self.groupBox_31.setSizePolicy(sizePolicy)
        self.groupBox_31.setMinimumSize(QtCore.QSize(570, 230))
        self.groupBox_31.setMaximumSize(QtCore.QSize(690, 250))
        font = QtGui.QFont()
        font.setFamily("Segoe UI Semibold")
        font.setPointSize(10)
        self.groupBox_31.setFont(font)
        self.groupBox_31.setStyleSheet("QGroupBox {\n"
"  \n"
"  border: 1px solid rgb(220, 220, 220);\n"
"  border-radius: 10px;\n"
"  padding: 4px;\n"
"  margin-top: 16px;\n"
"}\n"
"\n"
"QGroupBox::title {\n"
"  subcontrol-origin: margin;\n"
"  subcontrol-position: top left;\n"
"  left: 10px;\n"
"  padding-left: 3px;\n"
"  padding-right: 5px;\n"
"  padding-top: 8px;\n"
"  padding-bottom: 16px;\n"
"}\n"
"\n"
"QGroupBox::indicator {\n"
"  margin-left: 2px;\n"
"  height: 16px;\n"
"  width: 16px;\n"
"}")
        self.groupBox_31.setObjectName("groupBox_31")
        self.horizontalLayout_36 = QtWidgets.QHBoxLayout(self.groupBox_31)
        self.horizontalLayout_36.setObjectName("horizontalLayout_36")
        self.gridLayout_27 = QtWidgets.QGridLayout()
        self.gridLayout_27.setContentsMargins(3, 4, 0, 4)
        self.gridLayout_27.setHorizontalSpacing(3)
        self.gridLayout_27.setObjectName("gridLayout_27")
        self.label_102 = QtWidgets.QLabel(self.groupBox_31)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_102.sizePolicy().hasHeightForWidth())
        self.label_102.setSizePolicy(sizePolicy)
        self.label_102.setMinimumSize(QtCore.QSize(110, 0))
        self.label_102.setMaximumSize(QtCore.QSize(120, 16777215))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.label_102.setFont(font)
        self.label_102.setStyleSheet("QLabel{\n"
"margin-right: 20px;\n"
"}")
        self.label_102.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignVCenter)
        self.label_102.setObjectName("label_102")
        self.gridLayout_27.addWidget(self.label_102, 2, 0, 1, 1)
        self.label_133 = QtWidgets.QLabel(self.groupBox_31)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_133.sizePolicy().hasHeightForWidth())
        self.label_133.setSizePolicy(sizePolicy)
        self.label_133.setMinimumSize(QtCore.QSize(110, 0))
        self.label_133.setMaximumSize(QtCore.QSize(120, 16777215))
        font = QtGui.QFont()
        font.setPointSize(10)
        font.setWeight(50)
        font.setBold(False)
        self.label_133.setFont(font)
        self.label_133.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignVCenter)
        self.label_133.setObjectName("label_133")
        self.gridLayout_27.addWidget(self.label_133, 1, 0, 1, 1)
        self.conservation_PDB_ID_lineEdit = QtWidgets.QLineEdit(self.groupBox_31)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.conservation_PDB_ID_lineEdit.sizePolicy().hasHeightForWidth())
        self.conservation_PDB_ID_lineEdit.setSizePolicy(sizePolicy)
        self.conservation_PDB_ID_lineEdit.setMinimumSize(QtCore.QSize(100, 33))
        self.conservation_PDB_ID_lineEdit.setMaximumSize(QtCore.QSize(100, 33))
        self.conservation_PDB_ID_lineEdit.setStyleSheet("QLineEdit {\n"
"    background-color: rgb(27, 29, 35);\n"
"    border-radius: 5px;\n"
"    padding-left: 8px;\n"
"}\n"
"QLineEdit:hover {\n"
"    border: 2px solid rgb(64, 71, 88);\n"
"}\n"
"QLineEdit:focus {\n"
"    border: 2px solid rgb(91, 101, 124);\n"
"}")
        self.conservation_PDB_ID_lineEdit.setObjectName("conservation_PDB_ID_lineEdit")
        self.gridLayout_27.addWidget(self.conservation_PDB_ID_lineEdit, 0, 2, 1, 1)
        self.get_conserv_score_pushButton = QtWidgets.QPushButton(self.groupBox_31)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.get_conserv_score_pushButton.sizePolicy().hasHeightForWidth())
        self.get_conserv_score_pushButton.setSizePolicy(sizePolicy)
        self.get_conserv_score_pushButton.setMinimumSize(QtCore.QSize(120, 38))
        self.get_conserv_score_pushButton.setMaximumSize(QtCore.QSize(120, 38))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.get_conserv_score_pushButton.setFont(font)
        self.get_conserv_score_pushButton.setStyleSheet("QPushButton {\n"
"    border: 2px solid rgb(52, 59, 72);\n"
"    border-radius: 5px;    \n"
"    background-color: rgb(157, 90, 198);\n"
"    margin-top:1px;\n"
"    margin-bottom: 1px;\n"
"\n"
"}\n"
"QPushButton:hover {\n"
"    background-color: rgb(220, 5, 100);\n"
"    border: 2px solid rgb(61, 70, 86);\n"
"}\n"
"QPushButton:pressed {    \n"
"    background-color: rgb(190, 5, 100);\n"
"    border: 2px solid rgb(43, 50, 61);\n"
"}")
        self.get_conserv_score_pushButton.setObjectName("get_conserv_score_pushButton")
        self.gridLayout_27.addWidget(self.get_conserv_score_pushButton, 3, 0, 1, 1)
        self.conserv_score_doubleSpinBox = QtWidgets.QDoubleSpinBox(self.groupBox_31)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.conserv_score_doubleSpinBox.sizePolicy().hasHeightForWidth())
        self.conserv_score_doubleSpinBox.setSizePolicy(sizePolicy)
        self.conserv_score_doubleSpinBox.setMinimumSize(QtCore.QSize(100, 33))
        self.conserv_score_doubleSpinBox.setMaximumSize(QtCore.QSize(100, 33))
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(10)
        self.conserv_score_doubleSpinBox.setFont(font)
        self.conserv_score_doubleSpinBox.setLayoutDirection(QtCore.Qt.RightToLeft)
        self.conserv_score_doubleSpinBox.setStyleSheet("QDoubleSpinBox{\n"
"    background-color: rgb(27, 29, 35);\n"
"    border-radius: 5px;\n"
"    border: 2px solid rgb(27, 29, 35);\n"
"    padding: 5px;\n"
"    padding-left: 5px;\n"
"}\n"
"QDoubleSpinBox:hover{\n"
"    border: 2px solid rgb(64, 71, 88);\n"
"}\n"
"QDoubleSpinBox QAbstractItemView {\n"
"    color: rgb(85, 170, 255);    \n"
"    background-color: rgb(27, 29, 35);\n"
"    padding: 5px;\n"
"    selection-background-color: rgb(39, 44, 54);\n"
"}")
        self.conserv_score_doubleSpinBox.setDecimals(1)
        self.conserv_score_doubleSpinBox.setMinimum(1.0)
        self.conserv_score_doubleSpinBox.setMaximum(9.9)
        self.conserv_score_doubleSpinBox.setSingleStep(0.1)
        self.conserv_score_doubleSpinBox.setObjectName("conserv_score_doubleSpinBox")
        self.gridLayout_27.addWidget(self.conserv_score_doubleSpinBox, 2, 2, 1, 1)
        self.label_103 = QtWidgets.QLabel(self.groupBox_31)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_103.sizePolicy().hasHeightForWidth())
        self.label_103.setSizePolicy(sizePolicy)
        self.label_103.setMinimumSize(QtCore.QSize(110, 0))
        self.label_103.setMaximumSize(QtCore.QSize(120, 16777215))
        font = QtGui.QFont()
        font.setPointSize(10)
        font.setWeight(50)
        font.setBold(False)
        self.label_103.setFont(font)
        self.label_103.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignVCenter)
        self.label_103.setObjectName("label_103")
        self.gridLayout_27.addWidget(self.label_103, 0, 0, 1, 1)
        self.conservation_pdb_chain_id_lineedit = QtWidgets.QLineEdit(self.groupBox_31)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.conservation_pdb_chain_id_lineedit.sizePolicy().hasHeightForWidth())
        self.conservation_pdb_chain_id_lineedit.setSizePolicy(sizePolicy)
        self.conservation_pdb_chain_id_lineedit.setMinimumSize(QtCore.QSize(100, 33))
        self.conservation_pdb_chain_id_lineedit.setMaximumSize(QtCore.QSize(100, 33))
        self.conservation_pdb_chain_id_lineedit.setStyleSheet("QLineEdit {\n"
"    background-color: rgb(27, 29, 35);\n"
"    border-radius: 5px;\n"
"    padding-left: 8px;\n"
"\n"
"}\n"
"QLineEdit:hover {\n"
"    border: 2px solid rgb(64, 71, 88);\n"
"}\n"
"QLineEdit:focus {\n"
"    border: 2px solid rgb(91, 101, 124);\n"
"}")
        self.conservation_pdb_chain_id_lineedit.setObjectName("conservation_pdb_chain_id_lineedit")
        self.gridLayout_27.addWidget(self.conservation_pdb_chain_id_lineedit, 1, 2, 1, 1)
        self.verticalLayout_36 = QtWidgets.QVBoxLayout()
        self.verticalLayout_36.setSpacing(5)
        self.verticalLayout_36.setContentsMargins(3, 0, 3, 4)
        self.verticalLayout_36.setObjectName("verticalLayout_36")
        self.label_104 = QtWidgets.QLabel(self.groupBox_31)
        self.label_104.setEnabled(True)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_104.sizePolicy().hasHeightForWidth())
        self.label_104.setSizePolicy(sizePolicy)
        self.label_104.setMinimumSize(QtCore.QSize(300, 0))
        self.label_104.setMaximumSize(QtCore.QSize(16777215, 16777215))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.label_104.setFont(font)
        self.label_104.setStyleSheet("QLabel {\n"
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
        self.label_104.setAlignment(QtCore.Qt.AlignCenter)
        self.label_104.setObjectName("label_104")
        self.verticalLayout_36.addWidget(self.label_104)
        self.residues_conservation_tableWidget = QtWidgets.QTableWidget(self.groupBox_31)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.residues_conservation_tableWidget.sizePolicy().hasHeightForWidth())
        self.residues_conservation_tableWidget.setSizePolicy(sizePolicy)
        self.residues_conservation_tableWidget.setMinimumSize(QtCore.QSize(300, 150))
        self.residues_conservation_tableWidget.setMaximumSize(QtCore.QSize(16777215, 16777215))
        self.residues_conservation_tableWidget.setStyleSheet("\n"
"QTableView\n"
"{\n"
"    border-radius: 2px;\n"
"    border: 2px solid rgb(27, 29, 35);\n"
"    gridline-color: #6c6c6c;\n"
"    background-color:  rgb(27, 29, 35);\n"
"\n"
"    padding: 1px 1px 1px 1px;\n"
"    selection-background-color: rgb(127, 5, 64);\n"
"    selection-color: rgb(255, 255, 255);\n"
"    font-size: 10px;\n"
"}\n"
"\n"
"QWidget:hover  {\n"
"\n"
"    border: 2px solid rgb(64, 71, 88);\n"
"    selection-color: rgb(127, 5, 64);\n"
"}\n"
"\n"
"QTableView, QHeaderView\n"
"{\n"
"    border-radius: 0px;\n"
"}\n"
"\n"
"\n"
"QTableView::item:pressed, QListView::item:pressed, QTreeView::item:pressed  {\n"
"    background: #78879b;\n"
"    color: #FFFFFF;\n"
"}\n"
"\n"
"QTableView::item:selected:active, QTreeView::item:selected:active, QListView::item:selected:active  {\n"
"    background: #3d8ec9;\n"
"    color: #FFFFFF;\n"
"}\n"
"\n"
"\n"
"QHeaderView\n"
"{\n"
"    border: 1px transparent;\n"
"    border-radius: 2px;\n"
"    margin: 0px;\n"
"    padding: 3px;\n"
"}\n"
"\n"
"QHeaderView::section  {\n"
"    background-color: #3A3939;\n"
"    color: silver;\n"
"    padding: 2px;\n"
"    border: 1px solid #6c6c6c;\n"
"    border-radius: 2px;\n"
"    text-align: center;\n"
"}\n"
"\n"
"QHeaderView::section::vertical::first, QHeaderView::section::vertical::only-one\n"
"{\n"
"    border-top: 0px solid #6c6c6c;\n"
"}\n"
"\n"
"QHeaderView::section::vertical\n"
"{\n"
"    border-top: transparent;\n"
"}\n"
"\n"
"QHeaderView::section::horizontal::first, QHeaderView::section::horizontal::only-one\n"
"{\n"
"    border-left: 1px solid #6c6c6c;\n"
"}\n"
"\n"
"QHeaderView::section::horizontal\n"
"{\n"
"    border-left: transparent;\n"
"}\n"
"\n"
"\n"
"QHeaderView::section:checked\n"
" {\n"
"    color: white;\n"
"    background-color: #5A5959;\n"
" }\n"
"\n"
" /* style the sort indicator */\n"
"QHeaderView::down-arrow {\n"
"    image: url(:/16x16/icons/16x16/cil-arrow-circle-bottom.png);\n"
"}\n"
"\n"
"QHeaderView::up-arrow {\n"
"    image: url(:/16x16/icons/16x16/cil-arrow-circle-top.png);\n"
"}\n"
"\n"
"QTableCornerButton::section {\n"
"    background-color: #3A3939;\n"
"    border: 1px solid #3A3939;\n"
"    border-radius: 2px;\n"
"}\n"
"\n"
"\n"
"\n"
"QScrollBar:horizontal {\n"
"  border-color: black;\n"
"\n"
"  border-radius: 5px;\n"
"  background-color: black;\n"
"  height: 15px;\n"
"  margin: 0px 0px 1px 0px;\n"
" }\n"
"\n"
" QScrollBar::handle:horizontal {\n"
"    background-color:cyan;\n"
"    min-width: 15px;\n"
"    border-radius: 5px;\n"
" }\n"
"\n"
"\n"
"QScrollBar::add-line:horizontal {\n"
"    border-radius: 5px;\n"
"      background-color: rgb(58, 57, 57);\n"
"    height: 10px;\n"
"    subcontrol-position: right;\n"
"    subcontrol-origin: margin;\n"
" }\n"
"\n"
"\n"
" QScrollBar::sub-line:horizontal {\n"
"    border: 1px solid grey;\n"
"    background-color: rgb(241, 241, 241);\n"
"    height: 10px;\n"
"    subcontrol-position: left;\n"
"    subcontrol-origin: margin;\n"
" }\n"
"\n"
"\n"
"QScrollBar::add-page:horizontal, QScrollBar::sub-page:horizontal {\n"
"     background: none;\n"
"}\n"
"\n"
"")
        self.residues_conservation_tableWidget.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.residues_conservation_tableWidget.setWordWrap(True)
        self.residues_conservation_tableWidget.setObjectName("residues_conservation_tableWidget")
        self.residues_conservation_tableWidget.setColumnCount(2)
        self.residues_conservation_tableWidget.setRowCount(0)
        item = QtWidgets.QTableWidgetItem()
        self.residues_conservation_tableWidget.setHorizontalHeaderItem(0, item)
        item = QtWidgets.QTableWidgetItem()
        self.residues_conservation_tableWidget.setHorizontalHeaderItem(1, item)
        self.residues_conservation_tableWidget.horizontalHeader().setVisible(False)
        self.residues_conservation_tableWidget.horizontalHeader().setCascadingSectionResizes(False)
        self.residues_conservation_tableWidget.horizontalHeader().setDefaultSectionSize(150)
        self.residues_conservation_tableWidget.horizontalHeader().setHighlightSections(True)
        self.residues_conservation_tableWidget.horizontalHeader().setSortIndicatorShown(False)
        self.residues_conservation_tableWidget.horizontalHeader().setStretchLastSection(True)
        self.verticalLayout_36.addWidget(self.residues_conservation_tableWidget)
        self.gridLayout_27.addLayout(self.verticalLayout_36, 0, 4, 5, 1)
        self.use_conservation_checkBox = QtWidgets.QCheckBox(self.groupBox_31)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.use_conservation_checkBox.sizePolicy().hasHeightForWidth())
        self.use_conservation_checkBox.setSizePolicy(sizePolicy)
        self.use_conservation_checkBox.setMinimumSize(QtCore.QSize(115, 0))
        self.use_conservation_checkBox.setMaximumSize(QtCore.QSize(16777215, 16777215))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.use_conservation_checkBox.setFont(font)
        self.use_conservation_checkBox.setStyleSheet("QCheckBox {\n"
"    spacing: 5px;\n"
"    margin-left: 3px;\n"
"}\n"
"\n"
"QCheckBox::indicator {\n"
"    width: 15px;\n"
"    height: 15px;\n"
"}\n"
"\n"
"QCheckBox::indicator:unchecked {\n"
"border-image: url(:/16x16/icons/16x16/cil-circle.png);\n"
"\n"
"}\n"
"\n"
"\n"
"QCheckBox::indicator:checked {\n"
"    border-image: url(:/20x20/icons/20x20/cil-check.png);\n"
"}")
        self.use_conservation_checkBox.setObjectName("use_conservation_checkBox")
        self.gridLayout_27.addWidget(self.use_conservation_checkBox, 4, 0, 1, 2)
        self.horizontalLayout_36.addLayout(self.gridLayout_27)
        self.formLayout_2.setWidget(3, QtWidgets.QFormLayout.SpanningRole, self.groupBox_31)
        self.groupBox_30 = QtWidgets.QGroupBox(self.Perturb_Quick_3)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.groupBox_30.sizePolicy().hasHeightForWidth())
        self.groupBox_30.setSizePolicy(sizePolicy)
        self.groupBox_30.setMinimumSize(QtCore.QSize(610, 230))
        self.groupBox_30.setMaximumSize(QtCore.QSize(690, 230))
        font = QtGui.QFont()
        font.setFamily("Segoe UI Semibold")
        font.setPointSize(10)
        self.groupBox_30.setFont(font)
        self.groupBox_30.setStyleSheet("QGroupBox {\n"
"  \n"
"  border: 1px solid rgb(220, 220, 220);\n"
"  border-radius: 10px;\n"
"  padding: 4px;\n"
"  margin-top: 16px;\n"
"}\n"
"\n"
"QGroupBox::title {\n"
"  subcontrol-origin: margin;\n"
"  subcontrol-position: top left;\n"
"  left: 10px;\n"
"  padding-left: 3px;\n"
"  padding-right: 5px;\n"
"  padding-top: 8px;\n"
"  padding-bottom: 16px;\n"
"}\n"
"\n"
"QGroupBox::indicator {\n"
"  margin-left: 2px;\n"
"  height: 16px;\n"
"  width: 16px;\n"
"}")
        self.groupBox_30.setObjectName("groupBox_30")
        self.horizontalLayout_37 = QtWidgets.QHBoxLayout(self.groupBox_30)
        self.horizontalLayout_37.setObjectName("horizontalLayout_37")
        self.gridLayout_21 = QtWidgets.QGridLayout()
        self.gridLayout_21.setContentsMargins(3, 4, 0, 4)
        self.gridLayout_21.setHorizontalSpacing(3)
        self.gridLayout_21.setObjectName("gridLayout_21")
        self.label_98 = QtWidgets.QLabel(self.groupBox_30)
        self.label_98.setMinimumSize(QtCore.QSize(120, 0))
        self.label_98.setMaximumSize(QtCore.QSize(130, 16777215))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.label_98.setFont(font)
        self.label_98.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignVCenter)
        self.label_98.setObjectName("label_98")
        self.gridLayout_21.addWidget(self.label_98, 3, 0, 1, 1)
        self.horizontalLayout_29 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_29.setSpacing(5)
        self.horizontalLayout_29.setContentsMargins(2, -1, 1, -1)
        self.horizontalLayout_29.setObjectName("horizontalLayout_29")
        self.discard_residue_from_targets_pushButton = QtWidgets.QPushButton(self.groupBox_30)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.discard_residue_from_targets_pushButton.sizePolicy().hasHeightForWidth())
        self.discard_residue_from_targets_pushButton.setSizePolicy(sizePolicy)
        self.discard_residue_from_targets_pushButton.setMinimumSize(QtCore.QSize(50, 33))
        self.discard_residue_from_targets_pushButton.setMaximumSize(QtCore.QSize(65, 33))
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(9)
        self.discard_residue_from_targets_pushButton.setFont(font)
        self.discard_residue_from_targets_pushButton.setStyleSheet("QPushButton {\n"
"    border: 2px solid rgb(52, 59, 72);\n"
"    border-radius: 5px;    \n"
"    background-color: rgb(0, 180, 95);\n"
"}\n"
"QPushButton:hover {\n"
"    background-color: rgb(0, 152, 84);\n"
"    border: 2px solid rgb(61, 70, 86);\n"
"}\n"
"QPushButton:pressed {    \n"
"    background-color: rgb(220, 5, 100);\n"
"    border: 2px solid rgb(43, 50, 61);\n"
"}")
        self.discard_residue_from_targets_pushButton.setText("")
        self.discard_residue_from_targets_pushButton.setIcon(icon4)
        self.discard_residue_from_targets_pushButton.setObjectName("discard_residue_from_targets_pushButton")
        self.horizontalLayout_29.addWidget(self.discard_residue_from_targets_pushButton)
        self.add_residue_to_targets_pushButton = QtWidgets.QPushButton(self.groupBox_30)
        self.add_residue_to_targets_pushButton.setMinimumSize(QtCore.QSize(50, 33))
        self.add_residue_to_targets_pushButton.setMaximumSize(QtCore.QSize(60, 33))
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(9)
        self.add_residue_to_targets_pushButton.setFont(font)
        self.add_residue_to_targets_pushButton.setStyleSheet("QPushButton {\n"
"    border: 2px solid rgb(52, 59, 72);\n"
"    border-radius: 5px;    \n"
"    background-color: rgb(0, 180, 95);\n"
"}\n"
"QPushButton:hover {\n"
"    background-color: rgb(0, 152, 84);\n"
"    border: 2px solid rgb(61, 70, 86);\n"
"}\n"
"QPushButton:pressed {    \n"
"    background-color: rgb(220, 5, 100);\n"
"    border: 2px solid rgb(43, 50, 61);\n"
"}")
        self.add_residue_to_targets_pushButton.setText("")
        self.add_residue_to_targets_pushButton.setIcon(icon5)
        self.add_residue_to_targets_pushButton.setObjectName("add_residue_to_targets_pushButton")
        self.horizontalLayout_29.addWidget(self.add_residue_to_targets_pushButton)
        self.gridLayout_21.addLayout(self.horizontalLayout_29, 0, 2, 1, 1)
        self.node_threshold_spinBox = QtWidgets.QSpinBox(self.groupBox_30)
        self.node_threshold_spinBox.setEnabled(False)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.node_threshold_spinBox.sizePolicy().hasHeightForWidth())
        self.node_threshold_spinBox.setSizePolicy(sizePolicy)
        self.node_threshold_spinBox.setMinimumSize(QtCore.QSize(140, 33))
        self.node_threshold_spinBox.setMaximumSize(QtCore.QSize(16777215, 33))
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(10)
        self.node_threshold_spinBox.setFont(font)
        self.node_threshold_spinBox.setLayoutDirection(QtCore.Qt.RightToLeft)
        self.node_threshold_spinBox.setStyleSheet("QSpinBox{\n"
"    background-color: rgb(27, 29, 35);\n"
"    border-radius: 5px;\n"
"    border: 2px solid rgb(27, 29, 35);\n"
"    padding: 5px;\n"
"    padding-left: 5px;\n"
"}\n"
"QSpinBox:hover{\n"
"    border: 2px solid rgb(64, 71, 88);\n"
"}\n"
"QSpinBox QAbstractItemView {\n"
"    color: rgb(85, 170, 255);    \n"
"    background-color: rgb(27, 29, 35);\n"
"    padding: 5px;\n"
"    selection-background-color: rgb(39, 44, 54);\n"
"}\n"
"\n"
"QSpinBox:disabled {\n"
"  background-color: #19232D;\n"
"  color: #787878;\n"
"  selection-background-color: #14506E;\n"
"  selection-color: #787878;\n"
"}")
        self.node_threshold_spinBox.setMinimum(2)
        self.node_threshold_spinBox.setObjectName("node_threshold_spinBox")
        self.gridLayout_21.addWidget(self.node_threshold_spinBox, 3, 1, 1, 1)
        self.target_res_comboBox = QtWidgets.QComboBox(self.groupBox_30)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.target_res_comboBox.sizePolicy().hasHeightForWidth())
        self.target_res_comboBox.setSizePolicy(sizePolicy)
        self.target_res_comboBox.setMinimumSize(QtCore.QSize(110, 33))
        self.target_res_comboBox.setMaximumSize(QtCore.QSize(16777215, 33))
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(10)
        self.target_res_comboBox.setFont(font)
        self.target_res_comboBox.setStyleSheet("QComboBox {\n"
"    background-color: rgb(27, 29, 35);\n"
"    border-radius: 5px;\n"
"    border: 2px solid rgb(27, 29, 35);\n"
"    padding: 1px 18px 1px 3px;\n"
"    padding-left: 10px;\n"
"    min-width: 6em;\n"
"}\n"
"QComboBox:hover{\n"
"    border: 2px solid rgb(64, 71, 88);\n"
"}\n"
"\n"
"QComboBox QAbstractItemView{\n"
"    \n"
"    selection-background-color: rgb(200, 5, 100);\n"
"    \n"
"\n"
"}\n"
"\n"
"QComboBox:editable {\n"
"\n"
"    background: rgb(0, 180, 95);\n"
"}\n"
"\n"
"QComboBox:on { /* shift the text when the popup opens */\n"
"    padding-top: 3px;\n"
"    padding-left: 4px;\n"
"    text-align: center;\n"
"}\n"
"\n"
"\n"
"\n"
"QComboBox:on { /* shift the text when the popup opens */\n"
"    padding-top: 2px;\n"
"    padding-left: 3px;\n"
"}\n"
"\n"
"QComboBox::drop-down {\n"
"    subcontrol-origin: padding;\n"
"    subcontrol-position: top right;\n"
"    width: 20px;\n"
"\n"
"    border-left-width: 2px;\n"
"    border-left-color: darkgray;\n"
"    border-left-style: solid; /* just a single line */\n"
"    border-top-right-radius: 3px; /* same radius as the QComboBox */\n"
"    border-bottom-right-radius: 3px;\n"
"}\n"
"\n"
"QComboBox::down-arrow {\n"
"image: url(:/16x16/icons/16x16/cil-arrow-bottom.png);\n"
"}\n"
"\n"
"QComboBox::down-arrow:on { /* shift the arrow when popup is open */\n"
"    top: 1px;\n"
"    left: 1px;\n"
"}\n"
"\n"
"")
        self.target_res_comboBox.setObjectName("target_res_comboBox")
        self.gridLayout_21.addWidget(self.target_res_comboBox, 0, 1, 1, 1)
        self.network_cutoff_spinBox = QtWidgets.QSpinBox(self.groupBox_30)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.network_cutoff_spinBox.sizePolicy().hasHeightForWidth())
        self.network_cutoff_spinBox.setSizePolicy(sizePolicy)
        self.network_cutoff_spinBox.setMinimumSize(QtCore.QSize(0, 33))
        self.network_cutoff_spinBox.setMaximumSize(QtCore.QSize(16777215, 33))
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(10)
        self.network_cutoff_spinBox.setFont(font)
        self.network_cutoff_spinBox.setLayoutDirection(QtCore.Qt.RightToLeft)
        self.network_cutoff_spinBox.setStyleSheet("QSpinBox{\n"
"    background-color: rgb(27, 29, 35);\n"
"    border-radius: 5px;\n"
"    border: 2px solid rgb(27, 29, 35);\n"
"    padding: 5px;\n"
"    padding-left: 5px;\n"
"}\n"
"QSpinBox:hover{\n"
"    border: 2px solid rgb(64, 71, 88);\n"
"}\n"
"QSpinBox QAbstractItemView {\n"
"    color: rgb(85, 170, 255);    \n"
"    background-color: rgb(27, 29, 35);\n"
"    padding: 5px;\n"
"    selection-background-color: rgb(39, 44, 54);\n"
"}")
        self.network_cutoff_spinBox.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignVCenter)
        self.network_cutoff_spinBox.setButtonSymbols(QtWidgets.QAbstractSpinBox.UpDownArrows)
        self.network_cutoff_spinBox.setSpecialValueText("")
        self.network_cutoff_spinBox.setAccelerated(False)
        self.network_cutoff_spinBox.setCorrectionMode(QtWidgets.QAbstractSpinBox.CorrectToNearestValue)
        self.network_cutoff_spinBox.setMinimum(1)
        self.network_cutoff_spinBox.setMaximum(100)
        self.network_cutoff_spinBox.setProperty("value", 7)
        self.network_cutoff_spinBox.setObjectName("network_cutoff_spinBox")
        self.gridLayout_21.addWidget(self.network_cutoff_spinBox, 2, 1, 1, 1)
        self.label_99 = QtWidgets.QLabel(self.groupBox_30)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_99.sizePolicy().hasHeightForWidth())
        self.label_99.setSizePolicy(sizePolicy)
        self.label_99.setMinimumSize(QtCore.QSize(120, 0))
        self.label_99.setMaximumSize(QtCore.QSize(130, 16777215))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.label_99.setFont(font)
        self.label_99.setStyleSheet("QLabel{\n"
"margin-right: 20px;\n"
"}")
        self.label_99.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignVCenter)
        self.label_99.setObjectName("label_99")
        self.gridLayout_21.addWidget(self.label_99, 2, 0, 1, 1)
        self.node_threshold_checkBox = QtWidgets.QCheckBox(self.groupBox_30)
        self.node_threshold_checkBox.setMaximumSize(QtCore.QSize(115, 16777215))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.node_threshold_checkBox.setFont(font)
        self.node_threshold_checkBox.setStyleSheet("QCheckBox {\n"
"    spacing: 5px;\n"
"    margin-left: 15px;\n"
"}\n"
"\n"
"QCheckBox::indicator {\n"
"    width: 15px;\n"
"    height: 15px;\n"
"}\n"
"\n"
"QCheckBox::indicator:unchecked {\n"
"border-image: url(:/16x16/icons/16x16/cil-circle.png);\n"
"\n"
"}\n"
"\n"
"\n"
"QCheckBox::indicator:checked {\n"
"    border-image: url(:/20x20/icons/20x20/cil-check.png);\n"
"}")
        self.node_threshold_checkBox.setChecked(True)
        self.node_threshold_checkBox.setObjectName("node_threshold_checkBox")
        self.gridLayout_21.addWidget(self.node_threshold_checkBox, 3, 2, 1, 1)
        self.label_105 = QtWidgets.QLabel(self.groupBox_30)
        self.label_105.setMinimumSize(QtCore.QSize(120, 0))
        self.label_105.setMaximumSize(QtCore.QSize(130, 16777215))
        font = QtGui.QFont()
        font.setPointSize(10)
        font.setWeight(50)
        font.setBold(False)
        self.label_105.setFont(font)
        self.label_105.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignVCenter)
        self.label_105.setObjectName("label_105")
        self.gridLayout_21.addWidget(self.label_105, 1, 0, 1, 1)
        self.source_res_comboBox = QtWidgets.QComboBox(self.groupBox_30)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.source_res_comboBox.sizePolicy().hasHeightForWidth())
        self.source_res_comboBox.setSizePolicy(sizePolicy)
        self.source_res_comboBox.setMinimumSize(QtCore.QSize(110, 33))
        self.source_res_comboBox.setMaximumSize(QtCore.QSize(16777215, 33))
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(10)
        self.source_res_comboBox.setFont(font)
        self.source_res_comboBox.setStyleSheet("QComboBox {\n"
"    background-color: rgb(27, 29, 35);\n"
"    border-radius: 5px;\n"
"    border: 2px solid rgb(27, 29, 35);\n"
"    padding: 1px 18px 1px 3px;\n"
"    padding-left: 10px;\n"
"    min-width: 6em;\n"
"}\n"
"QComboBox:hover{\n"
"    border: 2px solid rgb(64, 71, 88);\n"
"}\n"
"\n"
"QComboBox QAbstractItemView{\n"
"    \n"
"    selection-background-color: rgb(200, 5, 100);\n"
"    \n"
"\n"
"}\n"
"\n"
"QComboBox:editable {\n"
"\n"
"    background: rgb(0, 180, 95);\n"
"}\n"
"\n"
"QComboBox:on { /* shift the text when the popup opens */\n"
"    padding-top: 3px;\n"
"    padding-left: 4px;\n"
"    text-align: center;\n"
"}\n"
"\n"
"\n"
"\n"
"QComboBox:on { /* shift the text when the popup opens */\n"
"    padding-top: 2px;\n"
"    padding-left: 3px;\n"
"}\n"
"\n"
"QComboBox::drop-down {\n"
"    subcontrol-origin: padding;\n"
"    subcontrol-position: top right;\n"
"    width: 20px;\n"
"\n"
"    border-left-width: 2px;\n"
"    border-left-color: darkgray;\n"
"    border-left-style: solid; /* just a single line */\n"
"    border-top-right-radius: 3px; /* same radius as the QComboBox */\n"
"    border-bottom-right-radius: 3px;\n"
"}\n"
"\n"
"QComboBox::down-arrow {\n"
"image: url(:/16x16/icons/16x16/cil-arrow-bottom.png);\n"
"}\n"
"\n"
"QComboBox::down-arrow:on { /* shift the arrow when popup is open */\n"
"    top: 1px;\n"
"    left: 1px;\n"
"}\n"
"\n"
"")
        self.source_res_comboBox.setObjectName("source_res_comboBox")
        self.gridLayout_21.addWidget(self.source_res_comboBox, 1, 1, 1, 1)
        self.label_97 = QtWidgets.QLabel(self.groupBox_30)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_97.sizePolicy().hasHeightForWidth())
        self.label_97.setSizePolicy(sizePolicy)
        self.label_97.setMinimumSize(QtCore.QSize(120, 0))
        self.label_97.setMaximumSize(QtCore.QSize(130, 16777215))
        font = QtGui.QFont()
        font.setPointSize(10)
        font.setWeight(50)
        font.setBold(False)
        self.label_97.setFont(font)
        self.label_97.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignVCenter)
        self.label_97.setObjectName("label_97")
        self.gridLayout_21.addWidget(self.label_97, 0, 0, 1, 1)
        self.verticalLayout_33 = QtWidgets.QVBoxLayout()
        self.verticalLayout_33.setSpacing(5)
        self.verticalLayout_33.setContentsMargins(3, 4, 3, 4)
        self.verticalLayout_33.setObjectName("verticalLayout_33")
        self.label_101 = QtWidgets.QLabel(self.groupBox_30)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_101.sizePolicy().hasHeightForWidth())
        self.label_101.setSizePolicy(sizePolicy)
        self.label_101.setMaximumSize(QtCore.QSize(215, 16777215))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.label_101.setFont(font)
        self.label_101.setStyleSheet("QLabel {\n"
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
        self.label_101.setAlignment(QtCore.Qt.AlignCenter)
        self.label_101.setObjectName("label_101")
        self.verticalLayout_33.addWidget(self.label_101)
        self.selected_target_residues_listWidget = QtWidgets.QListWidget(self.groupBox_30)
        self.selected_target_residues_listWidget.setMinimumSize(QtCore.QSize(0, 150))
        self.selected_target_residues_listWidget.setMaximumSize(QtCore.QSize(215, 16777215))
        self.selected_target_residues_listWidget.setStyleSheet("QListWidget {\n"
"    background-color: rgb(27, 29, 35);\n"
"    border-radius: 5px;\n"
"    border: 2px solid rgb(27, 29, 35);\n"
"    padding: 1px 18px 1px 3px;\n"
"    padding-left: 10px;\n"
"    selection-background-color: rgb(127, 5, 64);\n"
"    selection-color: rgb(255, 255, 255);\n"
"    font-size: 13px;\n"
"}\n"
"\n"
"QListWidget:hover{\n"
"    border: 2px solid rgb(64, 71, 88);\n"
"    selection-color: rgb(127, 5, 64);\n"
"\n"
"}\n"
"\n"
" QScrollBar:vertical\n"
"        {\n"
"            background-color: #323232;\n"
"            width: 15px;\n"
"            margin: 15px 3px 15px 3px;\n"
"            border: 1px transparent #2A2929;\n"
"            margin-right: 5px;\n"
"        }\n"
"        \n"
"        QScrollBar::handle:vertical\n"
"        {\n"
"            background-color: rgb(255, 17, 100);\n"
"            min-height: 5px;\n"
"            border-radius: 3px;\n"
"            \n"
"        }\n"
"        \n"
"        QScrollBar::sub-line:vertical\n"
"        {\n"
"            border: none;\n"
"            background: none;\n"
"            color: none;\n"
"        }\n"
"        \n"
"        QScrollBar::add-line:vertical\n"
"        {\n"
"            border: none;\n"
"            background: none;\n"
"            color: none;\n"
"        }\n"
"        \n"
"        QScrollBar::add-page:vertical, QScrollBar::sub-page:vertical\n"
"        {\n"
"            background: none;\n"
"        }\n"
"\n"
"")
        self.selected_target_residues_listWidget.setObjectName("selected_target_residues_listWidget")
        self.verticalLayout_33.addWidget(self.selected_target_residues_listWidget)
        self.gridLayout_21.addLayout(self.verticalLayout_33, 0, 3, 4, 1)
        self.horizontalLayout_37.addLayout(self.gridLayout_21)
        self.formLayout_2.setWidget(2, QtWidgets.QFormLayout.SpanningRole, self.groupBox_30)
        self.groupBox_28 = QtWidgets.QGroupBox(self.Perturb_Quick_3)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.groupBox_28.sizePolicy().hasHeightForWidth())
        self.groupBox_28.setSizePolicy(sizePolicy)
        self.groupBox_28.setMinimumSize(QtCore.QSize(610, 130))
        self.groupBox_28.setMaximumSize(QtCore.QSize(690, 130))
        font = QtGui.QFont()
        font.setFamily("Segoe UI Semibold")
        font.setPointSize(10)
        self.groupBox_28.setFont(font)
        self.groupBox_28.setStyleSheet("QGroupBox {\n"
"  \n"
"  border: 1px solid rgb(220, 220, 220);\n"
"  border-radius: 10px;\n"
"  padding: 4px;\n"
"  margin-top: 16px;\n"
"}\n"
"\n"
"QGroupBox::title {\n"
"  subcontrol-origin: margin;\n"
"  subcontrol-position: top left;\n"
"  left: 10px;\n"
"  padding-left: 3px;\n"
"  padding-right: 5px;\n"
"  padding-top: 8px;\n"
"  padding-bottom: 16px;\n"
"}\n"
"\n"
"QGroupBox::indicator {\n"
"  margin-left: 2px;\n"
"  height: 16px;\n"
"  width: 16px;\n"
"}")
        self.groupBox_28.setObjectName("groupBox_28")
        self.verticalLayout_31 = QtWidgets.QVBoxLayout(self.groupBox_28)
        self.verticalLayout_31.setObjectName("verticalLayout_31")
        self.gridLayout_20 = QtWidgets.QGridLayout()
        self.gridLayout_20.setContentsMargins(3, 4, 0, 4)
        self.gridLayout_20.setHorizontalSpacing(6)
        self.gridLayout_20.setObjectName("gridLayout_20")
        self.label_94 = QtWidgets.QLabel(self.groupBox_28)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_94.sizePolicy().hasHeightForWidth())
        self.label_94.setSizePolicy(sizePolicy)
        self.label_94.setMinimumSize(QtCore.QSize(120, 0))
        self.label_94.setMaximumSize(QtCore.QSize(130, 16777215))
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(10)
        self.label_94.setFont(font)
        self.label_94.setObjectName("label_94")
        self.gridLayout_20.addWidget(self.label_94, 0, 0, 1, 1)
        self.label_95 = QtWidgets.QLabel(self.groupBox_28)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_95.sizePolicy().hasHeightForWidth())
        self.label_95.setSizePolicy(sizePolicy)
        self.label_95.setMinimumSize(QtCore.QSize(120, 0))
        self.label_95.setMaximumSize(QtCore.QSize(130, 16777215))
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(10)
        self.label_95.setFont(font)
        self.label_95.setObjectName("label_95")
        self.gridLayout_20.addWidget(self.label_95, 1, 0, 1, 1)
        self.output_directory_button = QtWidgets.QPushButton(self.groupBox_28)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.output_directory_button.sizePolicy().hasHeightForWidth())
        self.output_directory_button.setSizePolicy(sizePolicy)
        self.output_directory_button.setMinimumSize(QtCore.QSize(106, 38))
        self.output_directory_button.setMaximumSize(QtCore.QSize(106, 38))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.output_directory_button.setFont(font)
        self.output_directory_button.setStyleSheet("QPushButton {\n"
"    border: 2px solid rgb(52, 59, 72);\n"
"    border-radius: 5px;    \n"
"    background-color: rgb(0, 180, 95);\n"
"    margin-top:1px;\n"
"    margin-bottom: 1px;\n"
"\n"
"}\n"
"QPushButton:hover {\n"
"    background-color: rgb(0, 152, 84);\n"
"    border: 2px solid rgb(61, 70, 86);\n"
"}\n"
"QPushButton:pressed {    \n"
"    background-color: rgb(220, 5, 100);\n"
"    border: 2px solid rgb(43, 50, 61);\n"
"}")
        self.output_directory_button.setObjectName("output_directory_button")
        self.gridLayout_20.addWidget(self.output_directory_button, 0, 3, 1, 1)
        self.net_output_directory_lineedit = QtWidgets.QLineEdit(self.groupBox_28)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.net_output_directory_lineedit.sizePolicy().hasHeightForWidth())
        self.net_output_directory_lineedit.setSizePolicy(sizePolicy)
        self.net_output_directory_lineedit.setMinimumSize(QtCore.QSize(340, 33))
        self.net_output_directory_lineedit.setMaximumSize(QtCore.QSize(16777215, 33))
        self.net_output_directory_lineedit.setStyleSheet("QLineEdit {\n"
"    background-color: rgb(27, 29, 35);\n"
"    border-radius: 5px;\n"
"    padding-left: 8px;\n"
"}\n"
"QLineEdit:hover {\n"
"    border: 2px solid rgb(64, 71, 88);\n"
"}\n"
"QLineEdit:focus {\n"
"    border: 2px solid rgb(91, 101, 124);\n"
"}")
        self.net_output_directory_lineedit.setObjectName("net_output_directory_lineedit")
        self.gridLayout_20.addWidget(self.net_output_directory_lineedit, 0, 1, 1, 1)
        self.PPI_Network_name_lineedit = QtWidgets.QLineEdit(self.groupBox_28)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.PPI_Network_name_lineedit.sizePolicy().hasHeightForWidth())
        self.PPI_Network_name_lineedit.setSizePolicy(sizePolicy)
        self.PPI_Network_name_lineedit.setMinimumSize(QtCore.QSize(340, 33))
        self.PPI_Network_name_lineedit.setMaximumSize(QtCore.QSize(16777215, 33))
        self.PPI_Network_name_lineedit.setStyleSheet("QLineEdit {\n"
"    background-color: rgb(27, 29, 35);\n"
"    border-radius: 5px;\n"
"\n"
"}\n"
"QLineEdit:hover {\n"
"    border: 2px solid rgb(64, 71, 88);\n"
"}\n"
"QLineEdit:focus {\n"
"    border: 2px solid rgb(91, 101, 124);\n"
"}")
        self.PPI_Network_name_lineedit.setObjectName("PPI_Network_name_lineedit")
        self.gridLayout_20.addWidget(self.PPI_Network_name_lineedit, 1, 1, 1, 1)
        self.verticalLayout_31.addLayout(self.gridLayout_20)
        self.formLayout_2.setWidget(1, QtWidgets.QFormLayout.SpanningRole, self.groupBox_28)
        self.groupBox_27 = QtWidgets.QGroupBox(self.Perturb_Quick_3)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.groupBox_27.sizePolicy().hasHeightForWidth())
        self.groupBox_27.setSizePolicy(sizePolicy)
        self.groupBox_27.setMinimumSize(QtCore.QSize(610, 130))
        self.groupBox_27.setMaximumSize(QtCore.QSize(690, 130))
        font = QtGui.QFont()
        font.setFamily("Segoe UI Semibold")
        font.setPointSize(10)
        font.setUnderline(False)
        font.setStrikeOut(False)
        self.groupBox_27.setFont(font)
        self.groupBox_27.setStyleSheet("QGroupBox {\n"
"  \n"
"  border: 1px solid rgb(220, 220, 220);\n"
"  border-radius: 10px;\n"
"  padding: 4px;\n"
"    margin-top: 5px;\n"
"  \n"
"}\n"
"\n"
"QGroupBox::title {\n"
"  subcontrol-origin: margin;\n"
"  subcontrol-position: top left;\n"
"  left: 10px;\n"
"  padding-left: 3px;\n"
"  padding-right: 5px;\n"
"  padding-top: -5px;\n"
"  padding-bottom: 16px;\n"
"}\n"
"\n"
"QGroupBox::indicator {\n"
"  margin-left: 2px;\n"
"  height: 16px;\n"
"  width: 16px;\n"
"}")
        self.groupBox_27.setFlat(False)
        self.groupBox_27.setCheckable(False)
        self.groupBox_27.setObjectName("groupBox_27")
        self.verticalLayout_32 = QtWidgets.QVBoxLayout(self.groupBox_27)
        self.verticalLayout_32.setSpacing(6)
        self.verticalLayout_32.setContentsMargins(-1, 9, -1, -1)
        self.verticalLayout_32.setObjectName("verticalLayout_32")
        self.gridLayout_19 = QtWidgets.QGridLayout()
        self.gridLayout_19.setSizeConstraint(QtWidgets.QLayout.SetFixedSize)
        self.gridLayout_19.setContentsMargins(3, 4, 0, 4)
        self.gridLayout_19.setObjectName("gridLayout_19")
        self.label_93 = QtWidgets.QLabel(self.groupBox_27)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_93.sizePolicy().hasHeightForWidth())
        self.label_93.setSizePolicy(sizePolicy)
        self.label_93.setMinimumSize(QtCore.QSize(120, 0))
        self.label_93.setMaximumSize(QtCore.QSize(130, 16777215))
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(10)
        font.setWeight(50)
        font.setBold(False)
        self.label_93.setFont(font)
        self.label_93.setObjectName("label_93")
        self.gridLayout_19.addWidget(self.label_93, 0, 0, 1, 1)
        self.response_time_lineEdit = QtWidgets.QLineEdit(self.groupBox_27)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.response_time_lineEdit.sizePolicy().hasHeightForWidth())
        self.response_time_lineEdit.setSizePolicy(sizePolicy)
        self.response_time_lineEdit.setMinimumSize(QtCore.QSize(340, 33))
        self.response_time_lineEdit.setMaximumSize(QtCore.QSize(16777215, 33))
        self.response_time_lineEdit.setStyleSheet("QLineEdit {\n"
"    background-color: rgb(27, 29, 35);\n"
"    border-radius: 5px;\n"
"    padding-left: 8px;\n"
"}\n"
"QLineEdit:hover {\n"
"    border: 2px solid rgb(64, 71, 88);\n"
"}\n"
"QLineEdit:focus {\n"
"    border: 2px solid rgb(91, 101, 124);\n"
"}")
        self.response_time_lineEdit.setObjectName("response_time_lineEdit")
        self.gridLayout_19.addWidget(self.response_time_lineEdit, 0, 2, 1, 1)
        self.label_91 = QtWidgets.QLabel(self.groupBox_27)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_91.sizePolicy().hasHeightForWidth())
        self.label_91.setSizePolicy(sizePolicy)
        self.label_91.setMinimumSize(QtCore.QSize(120, 0))
        self.label_91.setMaximumSize(QtCore.QSize(130, 16777215))
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(10)
        self.label_91.setFont(font)
        self.label_91.setObjectName("label_91")
        self.gridLayout_19.addWidget(self.label_91, 1, 0, 1, 1)
        self.upload_boundForm_pdb_Button = QtWidgets.QPushButton(self.groupBox_27)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.upload_boundForm_pdb_Button.sizePolicy().hasHeightForWidth())
        self.upload_boundForm_pdb_Button.setSizePolicy(sizePolicy)
        self.upload_boundForm_pdb_Button.setMinimumSize(QtCore.QSize(106, 38))
        self.upload_boundForm_pdb_Button.setMaximumSize(QtCore.QSize(106, 38))
        font = QtGui.QFont()
        font.setPointSize(10)
        font.setWeight(50)
        font.setBold(False)
        self.upload_boundForm_pdb_Button.setFont(font)
        self.upload_boundForm_pdb_Button.setStyleSheet("QPushButton {\n"
"    border: 2px solid rgb(52, 59, 72);\n"
"    border-radius: 5px;    \n"
"    background-color: rgb(0, 180, 95);\n"
"    margin-top: 2px;\n"
"    margin-bottom: 2px;\n"
"}\n"
"QPushButton:hover {\n"
"    background-color: rgb(0, 152, 84);\n"
"    border: 2px solid rgb(61, 70, 86);\n"
"}\n"
"QPushButton:pressed {    \n"
"    background-color: rgb(220, 5, 100);\n"
"    border: 2px solid rgb(43, 50, 61);\n"
"}\n"
"\n"
"QPushButton:disabled {\n"
"background-color: rgb(93, 83, 91);\n"
"}")
        self.upload_boundForm_pdb_Button.setObjectName("upload_boundForm_pdb_Button")
        self.gridLayout_19.addWidget(self.upload_boundForm_pdb_Button, 1, 3, 1, 1)
        self.response_time_upload_Button = QtWidgets.QPushButton(self.groupBox_27)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.response_time_upload_Button.sizePolicy().hasHeightForWidth())
        self.response_time_upload_Button.setSizePolicy(sizePolicy)
        self.response_time_upload_Button.setMinimumSize(QtCore.QSize(0, 38))
        self.response_time_upload_Button.setMaximumSize(QtCore.QSize(16777215, 38))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.response_time_upload_Button.setFont(font)
        self.response_time_upload_Button.setStyleSheet("QPushButton {\n"
"    border: 2px solid rgb(52, 59, 72);\n"
"    border-radius: 5px;    \n"
"    background-color: rgb(0, 180, 95);\n"
"    margin-top: 2px;\n"
"    margin-bottom: 2px;\n"
"}\n"
"QPushButton:hover {\n"
"    background-color: rgb(0, 152, 84);\n"
"    border: 2px solid rgb(61, 70, 86);\n"
"}\n"
"QPushButton:pressed {    \n"
"    background-color: rgb(220, 5, 100);\n"
"    border: 2px solid rgb(43, 50, 61);\n"
"}")
        self.response_time_upload_Button.setObjectName("response_time_upload_Button")
        self.gridLayout_19.addWidget(self.response_time_upload_Button, 0, 3, 1, 1)
        self.boundForm_pdb_lineedit = QtWidgets.QLineEdit(self.groupBox_27)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.boundForm_pdb_lineedit.sizePolicy().hasHeightForWidth())
        self.boundForm_pdb_lineedit.setSizePolicy(sizePolicy)
        self.boundForm_pdb_lineedit.setMinimumSize(QtCore.QSize(340, 33))
        self.boundForm_pdb_lineedit.setMaximumSize(QtCore.QSize(16777215, 33))
        self.boundForm_pdb_lineedit.setStyleSheet("QLineEdit {\n"
"    background-color: rgb(27, 29, 35);\n"
"    border-radius: 5px;\n"
"    padding-left: 8px;\n"
"\n"
"}\n"
"QLineEdit:hover {\n"
"    border: 2px solid rgb(64, 71, 88);\n"
"}\n"
"QLineEdit:focus {\n"
"    border: 2px solid rgb(91, 101, 124);\n"
"}")
        self.boundForm_pdb_lineedit.setObjectName("boundForm_pdb_lineedit")
        self.gridLayout_19.addWidget(self.boundForm_pdb_lineedit, 1, 2, 1, 1)
        self.verticalLayout_32.addLayout(self.gridLayout_19)
        self.formLayout_2.setWidget(0, QtWidgets.QFormLayout.SpanningRole, self.groupBox_27)
        self.add_tabb = QtWidgets.QPushButton(self.Perturb_Quick_3)
        self.add_tabb.setObjectName("add_tabb")
        self.formLayout_2.setWidget(5, QtWidgets.QFormLayout.LabelRole, self.add_tabb)
        self.horizontalLayout_34.addLayout(self.formLayout_2)
        self.PyMOL_3D_network_Widget = QtWidgets.QFrame(self.Perturb_Quick_3)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(1)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.PyMOL_3D_network_Widget.sizePolicy().hasHeightForWidth())
        self.PyMOL_3D_network_Widget.setSizePolicy(sizePolicy)
        self.PyMOL_3D_network_Widget.setStyleSheet("QFrame\n"
"{\n"
"    border: 1px solid black;\n"
"    border-radius: 5px;\n"
"    border-top-color: rgb(157, 90, 198);\n"
"    border-left-color: rgb(157, 90, 198);\n"
"    border-bottom-color: rgb(157, 90, 198);\n"
"    border-right-color: rgb(157, 90, 198);\n"
"    margin-top: 5px\n"
" \n"
"}")
        self.PyMOL_3D_network_Widget.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.PyMOL_3D_network_Widget.setFrameShadow(QtWidgets.QFrame.Raised)
        self.PyMOL_3D_network_Widget.setObjectName("PyMOL_3D_network_Widget")
        self.horizontalLayout_34.addWidget(self.PyMOL_3D_network_Widget)
        self.analysis_TabWidget.addTab(self.Perturb_Quick_3, "")
        self.tab_5 = QtWidgets.QWidget()
        self.tab_5.setObjectName("tab_5")
        self.horizontalLayout_20 = QtWidgets.QHBoxLayout(self.tab_5)
        self.horizontalLayout_20.setObjectName("horizontalLayout_20")
        self.gridLayout_10 = QtWidgets.QGridLayout()
        self.gridLayout_10.setObjectName("gridLayout_10")
        self.label_28 = QtWidgets.QLabel(self.tab_5)
        self.label_28.setMinimumSize(QtCore.QSize(0, 22))
        self.label_28.setMaximumSize(QtCore.QSize(16777215, 22))
        self.label_28.setStyleSheet("QLabel {\n"
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
        self.label_28.setObjectName("label_28")
        self.gridLayout_10.addWidget(self.label_28, 2, 0, 1, 1)
        self.shortest_path_listWidget = QtWidgets.QListWidget(self.tab_5)
        self.shortest_path_listWidget.setMaximumSize(QtCore.QSize(500, 16777215))
        self.shortest_path_listWidget.setObjectName("shortest_path_listWidget")
        self.gridLayout_10.addWidget(self.shortest_path_listWidget, 1, 0, 1, 1)
        self.intersection_path_listWidget = QtWidgets.QListWidget(self.tab_5)
        self.intersection_path_listWidget.setMaximumSize(QtCore.QSize(500, 16777215))
        self.intersection_path_listWidget.setObjectName("intersection_path_listWidget")
        self.gridLayout_10.addWidget(self.intersection_path_listWidget, 3, 0, 1, 1)
        self.label_29 = QtWidgets.QLabel(self.tab_5)
        self.label_29.setMinimumSize(QtCore.QSize(0, 22))
        self.label_29.setMaximumSize(QtCore.QSize(16777215, 22))
        self.label_29.setStyleSheet("QLabel {\n"
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
        self.label_29.setObjectName("label_29")
        self.gridLayout_10.addWidget(self.label_29, 0, 0, 1, 1)
        self.pyMOL_3D_analysis_frame = QtWidgets.QFrame(self.tab_5)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.pyMOL_3D_analysis_frame.sizePolicy().hasHeightForWidth())
        self.pyMOL_3D_analysis_frame.setSizePolicy(sizePolicy)
        self.pyMOL_3D_analysis_frame.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.pyMOL_3D_analysis_frame.setFrameShadow(QtWidgets.QFrame.Raised)
        self.pyMOL_3D_analysis_frame.setObjectName("pyMOL_3D_analysis_frame")
        self.gridLayout_10.addWidget(self.pyMOL_3D_analysis_frame, 0, 1, 4, 1)
        self.horizontalLayout_20.addLayout(self.gridLayout_10)
        self.analysis_TabWidget.addTab(self.tab_5, "")
        self.horizontalLayout_33.addWidget(self.analysis_TabWidget)
        self.stackedWidget.addWidget(self.Analysis)
        self.page_settings = QtWidgets.QWidget()
        self.page_settings.setObjectName("page_settings")
        self.horizontalLayout_9 = QtWidgets.QHBoxLayout(self.page_settings)
        self.horizontalLayout_9.setObjectName("horizontalLayout_9")
        self.stackedWidget.addWidget(self.page_settings)
        self.page_contact = QtWidgets.QWidget()
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.page_contact.sizePolicy().hasHeightForWidth())
        self.page_contact.setSizePolicy(sizePolicy)
        self.page_contact.setObjectName("page_contact")
        self.verticalLayout_14 = QtWidgets.QVBoxLayout(self.page_contact)
        self.verticalLayout_14.setObjectName("verticalLayout_14")
        self.verticalLayout_13 = QtWidgets.QVBoxLayout()
        self.verticalLayout_13.setObjectName("verticalLayout_13")
        self.verticalLayout_12 = QtWidgets.QVBoxLayout()
        self.verticalLayout_12.setObjectName("verticalLayout_12")
        self.label_34 = QtWidgets.QLabel(self.page_contact)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_34.sizePolicy().hasHeightForWidth())
        self.label_34.setSizePolicy(sizePolicy)
        self.label_34.setMaximumSize(QtCore.QSize(16777215, 100))
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(40)
        self.label_34.setFont(font)
        self.label_34.setStyleSheet("")
        self.label_34.setAlignment(QtCore.Qt.AlignCenter)
        self.label_34.setObjectName("label_34")
        self.verticalLayout_12.addWidget(self.label_34)
        self.horizontalLayout_18 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_18.setObjectName("horizontalLayout_18")
        spacerItem2 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_18.addItem(spacerItem2)
        self.frame_2 = QtWidgets.QFrame(self.page_contact)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.frame_2.sizePolicy().hasHeightForWidth())
        self.frame_2.setSizePolicy(sizePolicy)
        self.frame_2.setMinimumSize(QtCore.QSize(275, 250))
        self.frame_2.setMaximumSize(QtCore.QSize(16777215, 16777215))
        self.frame_2.setStyleSheet("border-image: url(:/big_icons/icons/big_icons/style_icon_2.png);\n"
"background-color: rgb(255, 255, 255);\n"
"border: 2px;\n"
"border-radius:20px;\n"
"border-color: rgb(0, 255, 127);")
        self.frame_2.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_2.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_2.setObjectName("frame_2")
        self.horizontalLayout_18.addWidget(self.frame_2)
        spacerItem3 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_18.addItem(spacerItem3)
        self.verticalLayout_12.addLayout(self.horizontalLayout_18)
        self.label_33 = QtWidgets.QLabel(self.page_contact)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_33.sizePolicy().hasHeightForWidth())
        self.label_33.setSizePolicy(sizePolicy)
        self.label_33.setMinimumSize(QtCore.QSize(0, 25))
        self.label_33.setMaximumSize(QtCore.QSize(16777215, 100))
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(15)
        self.label_33.setFont(font)
        self.label_33.setAlignment(QtCore.Qt.AlignCenter)
        self.label_33.setObjectName("label_33")
        self.verticalLayout_12.addWidget(self.label_33)
        self.label_41 = QtWidgets.QLabel(self.page_contact)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_41.sizePolicy().hasHeightForWidth())
        self.label_41.setSizePolicy(sizePolicy)
        self.label_41.setMinimumSize(QtCore.QSize(0, 25))
        self.label_41.setMaximumSize(QtCore.QSize(16777215, 100))
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(15)
        self.label_41.setFont(font)
        self.label_41.setAlignment(QtCore.Qt.AlignCenter)
        self.label_41.setObjectName("label_41")
        self.verticalLayout_12.addWidget(self.label_41)
        self.textBrowser = QtWidgets.QTextBrowser(self.page_contact)
        self.textBrowser.setMinimumSize(QtCore.QSize(600, 0))
        self.textBrowser.setMaximumSize(QtCore.QSize(16777215, 300))
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        self.textBrowser.setFont(font)
        self.textBrowser.setStyleSheet("QTextBrowser {\n"
"    background-color: rgb(27, 29, 35);\n"
"    border-radius: 10px;\n"
"    padding-left: 20px;\n"
"    padding-right: 20px;\n"
"    margin-left: 200px;\n"
"    margin-right: 200px;\n"
"    margin-bottom: 10px;\n"
"    margin-top: 10px;\n"
"    \n"
"}\n"
"QTextBrowser:hover {\n"
"    border: 2px solid rgb(64, 71, 88);\n"
"}\n"
"QTextBrowser:focus {\n"
"    border: 2px solid rgb(91, 101, 124);\n"
"}\n"
"\n"
"\n"
"\n"
"        QScrollBar:vertical\n"
"        {\n"
"            background-color: #323232;\n"
"            width: 15px;\n"
"            margin: 15px 3px 15px 3px;\n"
"            border: 1px transparent #2A2929;\n"
"            margin-right: 5px;\n"
"        }\n"
"        \n"
"        QScrollBar::handle:vertical\n"
"        {\n"
"            background-color: rgb(255, 17, 100);\n"
"            min-height: 5px;\n"
"            border-radius: 3px;\n"
"            \n"
"        }\n"
"        \n"
"        QScrollBar::sub-line:vertical\n"
"        {\n"
"            border: none;\n"
"            background: none;\n"
"            color: none;\n"
"        }\n"
"        \n"
"        QScrollBar::add-line:vertical\n"
"        {\n"
"            border: none;\n"
"            background: none;\n"
"            color: none;\n"
"        }\n"
"        \n"
"        QScrollBar::add-page:vertical, QScrollBar::sub-page:vertical\n"
"        {\n"
"            background: none;\n"
"        }")
        self.textBrowser.setObjectName("textBrowser")
        self.verticalLayout_12.addWidget(self.textBrowser)
        self.horizontalLayout_17 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_17.setSpacing(20)
        self.horizontalLayout_17.setContentsMargins(200, 5, 6, 5)
        self.horizontalLayout_17.setObjectName("horizontalLayout_17")
        self.commandLinkButton_3 = QtWidgets.QCommandLinkButton(self.page_contact)
        self.commandLinkButton_3.setMinimumSize(QtCore.QSize(0, 100))
        self.commandLinkButton_3.setMaximumSize(QtCore.QSize(300, 30))
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(10)
        font.setWeight(75)
        font.setBold(True)
        self.commandLinkButton_3.setFont(font)
        self.commandLinkButton_3.setAutoFillBackground(False)
        self.commandLinkButton_3.setStyleSheet("QCommandLinkButton {        \n"
"    \n"
"    border-radius: 5px;\n"
"    color: rgb(255, 0, 127);\n"
"    \n"
"    text-align: center;\n"
"}\n"
"QCommandLinkButton:hover {    \n"
"    color: rgb(210, 210, 210);\n"
"    background-color: rgb(44, 49, 60);\n"
"    padding: 20px\n"
"}\n"
"QCommandLinkButton:pressed {    \n"
"    color: rgb(210, 210, 210);\n"
"    background-color: rgb(52, 58, 71);\n"
"}")
        icon19 = QtGui.QIcon()
        icon19.addPixmap(QtGui.QPixmap(":/24x24/icons/24x24/cil-link.png"), QtGui.QIcon.Normal, QtGui.QIcon.On)
        self.commandLinkButton_3.setIcon(icon19)
        self.commandLinkButton_3.setObjectName("commandLinkButton_3")
        self.horizontalLayout_17.addWidget(self.commandLinkButton_3)
        self.commandLinkButton_2 = QtWidgets.QCommandLinkButton(self.page_contact)
        self.commandLinkButton_2.setMaximumSize(QtCore.QSize(300, 30))
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(10)
        font.setWeight(75)
        font.setBold(True)
        self.commandLinkButton_2.setFont(font)
        self.commandLinkButton_2.setAutoFillBackground(False)
        self.commandLinkButton_2.setStyleSheet("QCommandLinkButton {        \n"
"    color: rgb(85, 170, 255);\n"
"    border-radius: 5px;\n"
"    padding: 5px;\n"
"}\n"
"QCommandLinkButton:hover {    \n"
"    color: rgb(210, 210, 210);\n"
"    background-color: rgb(44, 49, 60);\n"
"}\n"
"QCommandLinkButton:pressed {    \n"
"    color: rgb(210, 210, 210);\n"
"    background-color: rgb(52, 58, 71);\n"
"}")
        self.commandLinkButton_2.setIcon(icon19)
        self.commandLinkButton_2.setObjectName("commandLinkButton_2")
        self.horizontalLayout_17.addWidget(self.commandLinkButton_2)
        self.commandLinkButton = QtWidgets.QCommandLinkButton(self.page_contact)
        self.commandLinkButton.setMaximumSize(QtCore.QSize(300, 30))
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(10)
        font.setWeight(75)
        font.setBold(True)
        self.commandLinkButton.setFont(font)
        self.commandLinkButton.setAutoFillBackground(False)
        self.commandLinkButton.setStyleSheet("QCommandLinkButton {        \n"
"    color: rgb(85, 170, 255);\n"
"    border-radius: 5px;\n"
"    padding: 5px;\n"
"    \n"
"}\n"
"QCommandLinkButton:hover {    \n"
"    color: rgb(210, 210, 210);\n"
"    background-color: rgb(44, 49, 60);\n"
"}\n"
"QCommandLinkButton:pressed {    \n"
"    color: rgb(210, 210, 210);\n"
"    background-color: rgb(52, 58, 71);\n"
"}")
        self.commandLinkButton.setIcon(icon19)
        self.commandLinkButton.setObjectName("commandLinkButton")
        self.horizontalLayout_17.addWidget(self.commandLinkButton)
        self.verticalLayout_12.addLayout(self.horizontalLayout_17)
        self.verticalLayout_13.addLayout(self.verticalLayout_12)
        self.verticalLayout_14.addLayout(self.verticalLayout_13)
        self.stackedWidget.addWidget(self.page_contact)
        self.verticalLayout_9.addWidget(self.stackedWidget)
        self.verticalLayout_4.addWidget(self.frame_content)
        self.frame_grip = QtWidgets.QFrame(self.frame_content_right)
        self.frame_grip.setMinimumSize(QtCore.QSize(0, 25))
        self.frame_grip.setMaximumSize(QtCore.QSize(16777215, 25))
        self.frame_grip.setStyleSheet("background-color: rgb(33, 37, 43);")
        self.frame_grip.setFrameShape(QtWidgets.QFrame.NoFrame)
        self.frame_grip.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_grip.setObjectName("frame_grip")
        self.horizontalLayout_6 = QtWidgets.QHBoxLayout(self.frame_grip)
        self.horizontalLayout_6.setSpacing(0)
        self.horizontalLayout_6.setContentsMargins(0, 0, 2, 0)
        self.horizontalLayout_6.setObjectName("horizontalLayout_6")
        self.frame_label_bottom = QtWidgets.QFrame(self.frame_grip)
        self.frame_label_bottom.setFrameShape(QtWidgets.QFrame.NoFrame)
        self.frame_label_bottom.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_label_bottom.setObjectName("frame_label_bottom")
        self.horizontalLayout_7 = QtWidgets.QHBoxLayout(self.frame_label_bottom)
        self.horizontalLayout_7.setSpacing(0)
        self.horizontalLayout_7.setContentsMargins(10, 0, 10, 0)
        self.horizontalLayout_7.setObjectName("horizontalLayout_7")
        self.label_credits = QtWidgets.QLabel(self.frame_label_bottom)
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        self.label_credits.setFont(font)
        self.label_credits.setStyleSheet("color: rgb(98, 103, 111);")
        self.label_credits.setObjectName("label_credits")
        self.horizontalLayout_7.addWidget(self.label_credits)
        self.label_version = QtWidgets.QLabel(self.frame_label_bottom)
        self.label_version.setMaximumSize(QtCore.QSize(100, 16777215))
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        self.label_version.setFont(font)
        self.label_version.setStyleSheet("color: rgb(98, 103, 111);")
        self.label_version.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.label_version.setObjectName("label_version")
        self.horizontalLayout_7.addWidget(self.label_version)
        self.horizontalLayout_6.addWidget(self.frame_label_bottom)
        self.frame_size_grip = QtWidgets.QFrame(self.frame_grip)
        self.frame_size_grip.setMaximumSize(QtCore.QSize(20, 20))
        self.frame_size_grip.setStyleSheet("QSizeGrip {\n"
"    background-image: url(:/16x16/icons/16x16/cil-size-grip.png);\n"
"    background-position: center;\n"
"    background-repeat: no-reperat;\n"
"}")
        self.frame_size_grip.setFrameShape(QtWidgets.QFrame.NoFrame)
        self.frame_size_grip.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_size_grip.setObjectName("frame_size_grip")
        self.horizontalLayout_6.addWidget(self.frame_size_grip)
        self.verticalLayout_4.addWidget(self.frame_grip)
        self.horizontalLayout_2.addWidget(self.frame_content_right)
        self.verticalLayout.addWidget(self.frame_center)
        self.horizontalLayout.addWidget(self.frame_main)
        MainWindow.setCentralWidget(self.centralwidget)
        self.label_35.setBuddy(self.Number_CPU_spinBox)
        self.label_6.setBuddy(self.Number_CPU_spinBox)
        self.label_100.setBuddy(self.Number_CPU_spinBox)
        self.label_98.setBuddy(self.Number_CPU_spinBox)

        self.retranslateUi(MainWindow)
        self.stackedWidget.setCurrentIndex(2)
        self.Perturbation_TabWidget.setCurrentIndex(0)
        self.tabWidget.setCurrentIndex(0)
        self.analysis_TabWidget.setCurrentIndex(1)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        MainWindow.setWindowTitle(QtWidgets.QApplication.translate("MainWindow", "MainWindow", None, -1))
        self.label_title_bar_top.setText(QtWidgets.QApplication.translate("MainWindow", "MDPerTool", None, -1))
        self.btn_minimize.setToolTip(QtWidgets.QApplication.translate("MainWindow", "Minimize", None, -1))
        self.btn_maximize_restore.setToolTip(QtWidgets.QApplication.translate("MainWindow", "Maximize", None, -1))
        self.btn_close.setToolTip(QtWidgets.QApplication.translate("MainWindow", "Close", None, -1))
        self.label_top_info_1.setText(QtWidgets.QApplication.translate("MainWindow", "C:\\Program Files\\MDPerTool", None, -1))
        self.label_top_info_2.setText(QtWidgets.QApplication.translate("MainWindow", "| HOME", None, -1))
        self.label_user_icon.setText(QtWidgets.QApplication.translate("MainWindow", "HIO", None, -1))
        self.groupBox_14.setTitle(QtWidgets.QApplication.translate("MainWindow", "Inputs", None, -1))
        self.pdb_fix_clean.setText(QtWidgets.QApplication.translate("MainWindow", "PDB fix and clean", None, -1))
        self.label_32.setText(QtWidgets.QApplication.translate("MainWindow", " PDB File: ", None, -1))
        self.upload_pdb_Button.setText(QtWidgets.QApplication.translate("MainWindow", "Upload", None, -1))
        self.fetch_pdb_Button.setText(QtWidgets.QApplication.translate("MainWindow", "Fetch ", None, -1))
        self.label_2.setText(QtWidgets.QApplication.translate("MainWindow", "Or", None, -1))
        self.label.setText(QtWidgets.QApplication.translate("MainWindow", " PDB ID:", None, -1))
        self.groupBox_16.setTitle(QtWidgets.QApplication.translate("MainWindow", "Output", None, -1))
        self.label_31.setText(QtWidgets.QApplication.translate("MainWindow", " Folder: ", None, -1))
        self.Browse_Output_button.setText(QtWidgets.QApplication.translate("MainWindow", "Browse", None, -1))
        self.groupBox_17.setTitle(QtWidgets.QApplication.translate("MainWindow", "Long Simulation Parameters", None, -1))
        self.label_35.setText(QtWidgets.QApplication.translate("MainWindow", "Threading:", None, -1))
        self.All_CPU_checkBox.setText(QtWidgets.QApplication.translate("MainWindow", "All CPU", None, -1))
        self.long_simulation_time_unit.setItemText(0, QtWidgets.QApplication.translate("MainWindow", "nanosecond", None, -1))
        self.long_simulation_time_unit.setItemText(1, QtWidgets.QApplication.translate("MainWindow", "picosecond", None, -1))
        self.label_8.setText(QtWidgets.QApplication.translate("MainWindow", "Run Duration: ", None, -1))
        self.groupBox_15.setTitle(QtWidgets.QApplication.translate("MainWindow", "Short Simulation Parameters", None, -1))
        self.label_5.setText(QtWidgets.QApplication.translate("MainWindow", "Residue: ", None, -1))
        self.label_6.setText(QtWidgets.QApplication.translate("MainWindow", "Threading:", None, -1))
        self.checkBox_2.setText(QtWidgets.QApplication.translate("MainWindow", "All CPU", None, -1))
        self.perturb_time_unit_comboBox.setItemText(0, QtWidgets.QApplication.translate("MainWindow", "step", None, -1))
        self.label_4.setText(QtWidgets.QApplication.translate("MainWindow", "R Factor: ", None, -1))
        self.label_3.setText(QtWidgets.QApplication.translate("MainWindow", "Run Duration: ", None, -1))
        self.R_factor_lineEdit.setText(QtWidgets.QApplication.translate("MainWindow", "3,4", None, -1))
        self.label_7.setText(QtWidgets.QApplication.translate("MainWindow", "Selected Residue(s)", None, -1))
        self.quit_pushButton.setText(QtWidgets.QApplication.translate("MainWindow", "Quit", None, -1))
        self.stop_pushButton.setText(QtWidgets.QApplication.translate("MainWindow", "Stop", None, -1))
        self.pushButton_2.setText(QtWidgets.QApplication.translate("MainWindow", "Save Script", None, -1))
        self.Run.setText(QtWidgets.QApplication.translate("MainWindow", "Run", None, -1))
        self.Perturbation_TabWidget.setTabText(self.Perturbation_TabWidget.indexOf(self.Perturb_Quick), QtWidgets.QApplication.translate("MainWindow", "Quick", None, -1))
        self.groupBox_18.setTitle(QtWidgets.QApplication.translate("MainWindow", "Devices", None, -1))
        self.label_39.setText(QtWidgets.QApplication.translate("MainWindow", "Device ID: ", None, -1))
        self.platform_comboBox.setItemText(0, QtWidgets.QApplication.translate("MainWindow", "CUDA", None, -1))
        self.platform_comboBox.setItemText(1, QtWidgets.QApplication.translate("MainWindow", "OpenCL", None, -1))
        self.platform_comboBox.setItemText(2, QtWidgets.QApplication.translate("MainWindow", "Reference", None, -1))
        self.platform_comboBox.setItemText(3, QtWidgets.QApplication.translate("MainWindow", "CPU", None, -1))
        self.Device_Number_comboBox.setItemText(0, QtWidgets.QApplication.translate("MainWindow", "1", None, -1))
        self.label_38.setText(QtWidgets.QApplication.translate("MainWindow", "Platform: ", None, -1))
        self.Device_ID_checkBox.setText(QtWidgets.QApplication.translate("MainWindow", "Use Device ID", None, -1))
        self.groupBox_4.setTitle(QtWidgets.QApplication.translate("MainWindow", "Forcefield", None, -1))
        self.protein_forcefield_comboBox.setItemText(0, QtWidgets.QApplication.translate("MainWindow", "amber03", None, -1))
        self.protein_forcefield_comboBox.setItemText(1, QtWidgets.QApplication.translate("MainWindow", "amber10", None, -1))
        self.protein_forcefield_comboBox.setItemText(2, QtWidgets.QApplication.translate("MainWindow", "amber96", None, -1))
        self.protein_forcefield_comboBox.setItemText(3, QtWidgets.QApplication.translate("MainWindow", "amber99sb", None, -1))
        self.protein_forcefield_comboBox.setItemText(4, QtWidgets.QApplication.translate("MainWindow", "amber99sbildn", None, -1))
        self.protein_forcefield_comboBox.setItemText(5, QtWidgets.QApplication.translate("MainWindow", "charmm36", None, -1))
        self.label_10.setText(QtWidgets.QApplication.translate("MainWindow", "Protein: ", None, -1))
        self.label_11.setText(QtWidgets.QApplication.translate("MainWindow", "Water: ", None, -1))
        self.water_forcefield_comboBox.setItemText(0, QtWidgets.QApplication.translate("MainWindow", "tip3p", None, -1))
        self.water_forcefield_comboBox.setItemText(1, QtWidgets.QApplication.translate("MainWindow", "spce", None, -1))
        self.water_forcefield_comboBox.setItemText(2, QtWidgets.QApplication.translate("MainWindow", "tip5p", None, -1))
        self.water_forcefield_comboBox.setItemText(3, QtWidgets.QApplication.translate("MainWindow", "tip4pew", None, -1))
        self.label_42.setText(QtWidgets.QApplication.translate("MainWindow", "Padding:", None, -1))
        self.water_padding_lineEdit.setText(QtWidgets.QApplication.translate("MainWindow", "10 * angstrom", None, -1))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab), QtWidgets.QApplication.translate("MainWindow", "General", None, -1))
        self.groupBox_5.setTitle(QtWidgets.QApplication.translate("MainWindow", "Integrator Options", None, -1))
        self.groupBox_3.setTitle(QtWidgets.QApplication.translate("MainWindow", "GroupBox", None, -1))
        self.integrator_kind_comboBox.setItemText(0, QtWidgets.QApplication.translate("MainWindow", "Langevin", None, -1))
        self.integrator_kind_comboBox.setItemText(1, QtWidgets.QApplication.translate("MainWindow", "Brownian", None, -1))
        self.integrator_kind_comboBox.setItemText(2, QtWidgets.QApplication.translate("MainWindow", "Verlet", None, -1))
        self.integrator_time_step.setHtml(QtWidgets.QApplication.translate("MainWindow", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'Segoe UI\'; font-size:10pt; font-weight:400; font-style:normal;\">\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'MS Shell Dlg 2\'; font-size:8.25pt;\">2.0</span></p></body></html>", None, -1))
        self.label_12.setText(QtWidgets.QApplication.translate("MainWindow", "Timestep: ", None, -1))
        self.label_9.setText(QtWidgets.QApplication.translate("MainWindow", "Kind: ", None, -1))
        self.Additional_Integrators_checkBox.setText(QtWidgets.QApplication.translate("MainWindow", "Additional Integrators", None, -1))
        self.integrator_time_step_unit.setItemText(0, QtWidgets.QApplication.translate("MainWindow", "femtosecond", None, -1))
        self.Additional_Integrator_groupBox.setTitle(QtWidgets.QApplication.translate("MainWindow", "Additional Integrators", None, -1))
        self.label_13.setText(QtWidgets.QApplication.translate("MainWindow", "Friction: ", None, -1))
        self.label_14.setText(QtWidgets.QApplication.translate("MainWindow", "Temperature: ", None, -1))
        self.friction_textEdit.setHtml(QtWidgets.QApplication.translate("MainWindow", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'Segoe UI\'; font-size:10pt; font-weight:400; font-style:normal;\">\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'MS Shell Dlg 2\'; font-size:8.25pt;\">91.0/picosecond</span></p></body></html>", None, -1))
        self.temperature_textEdit.setHtml(QtWidgets.QApplication.translate("MainWindow", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'Segoe UI\'; font-size:10pt; font-weight:400; font-style:normal;\">\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'MS Shell Dlg 2\'; font-size:8.25pt;\">310.0*kelvin</span></p></body></html>", None, -1))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_2), QtWidgets.QApplication.translate("MainWindow", "Integrator", None, -1))
        self.groupBox_7.setTitle(QtWidgets.QApplication.translate("MainWindow", "System Options", None, -1))
        self.system_constraints_comboBox.setItemText(0, QtWidgets.QApplication.translate("MainWindow", "None", None, -1))
        self.system_constraints_comboBox.setItemText(1, QtWidgets.QApplication.translate("MainWindow", "HBonds", None, -1))
        self.system_constraints_comboBox.setItemText(2, QtWidgets.QApplication.translate("MainWindow", "HAngles", None, -1))
        self.system_constraints_comboBox.setItemText(3, QtWidgets.QApplication.translate("MainWindow", "AllBonds", None, -1))
        self.label_17.setText(QtWidgets.QApplication.translate("MainWindow", "Rigid water: ", None, -1))
        self.label_18.setText(QtWidgets.QApplication.translate("MainWindow", "Nonbounded cutoff: ", None, -1))
        self.label_15.setText(QtWidgets.QApplication.translate("MainWindow", "Nonbounded Method: ", None, -1))
        self.label_16.setText(QtWidgets.QApplication.translate("MainWindow", "Constraints: ", None, -1))
        self.nonbounded_CutOff_textEdit.setHtml(QtWidgets.QApplication.translate("MainWindow", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'Segoe UI\'; font-size:10pt; font-weight:400; font-style:normal;\">\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'MS Shell Dlg 2\'; font-size:8.25pt;\">1.2*nanometer</span></p></body></html>", None, -1))
        self.nonBounded_Method_comboBox.setItemText(0, QtWidgets.QApplication.translate("MainWindow", "PME", None, -1))
        self.nonBounded_Method_comboBox.setItemText(1, QtWidgets.QApplication.translate("MainWindow", "NoCutoff", None, -1))
        self.nonBounded_Method_comboBox.setItemText(2, QtWidgets.QApplication.translate("MainWindow", "CutoffNonPeriodic", None, -1))
        self.nonBounded_Method_comboBox.setItemText(3, QtWidgets.QApplication.translate("MainWindow", "CutoffPeriodic", None, -1))
        self.nonBounded_Method_comboBox.setItemText(4, QtWidgets.QApplication.translate("MainWindow", "Ewald", None, -1))
        self.switching_distance_textEdit.setHtml(QtWidgets.QApplication.translate("MainWindow", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'Segoe UI\'; font-size:10pt; font-weight:400; font-style:normal;\">\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'MS Shell Dlg 2\'; font-size:8.25pt;\">1*nanometer</span></p></body></html>", None, -1))
        self.label_45.setText(QtWidgets.QApplication.translate("MainWindow", "Switching Distance", None, -1))
        self.use_switching_checkBox.setText(QtWidgets.QApplication.translate("MainWindow", "Use Switching Distance", None, -1))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_3), QtWidgets.QApplication.translate("MainWindow", "System", None, -1))
        self.groupBox_8.setTitle(QtWidgets.QApplication.translate("MainWindow", "Simulation Options", None, -1))
        self.groupBox_10.setTitle(QtWidgets.QApplication.translate("MainWindow", "Minimization Options", None, -1))
        self.label_20.setText(QtWidgets.QApplication.translate("MainWindow", "Minimize: ", None, -1))
        self.label_21.setText(QtWidgets.QApplication.translate("MainWindow", "Max iterations: ", None, -1))
        self.Max_minimize_iter_textEdit.setHtml(QtWidgets.QApplication.translate("MainWindow", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'Segoe UI\'; font-size:10pt; font-weight:400; font-style:normal;\">\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'MS Shell Dlg 2\'; font-size:8.25pt;\">500</span></p></body></html>", None, -1))
        self.label_43.setText(QtWidgets.QApplication.translate("MainWindow", "Equilibrate:", None, -1))
        self.label_44.setText(QtWidgets.QApplication.translate("MainWindow", "Equilibrate steps: ", None, -1))
        self.Max_equilubrate_steps_textEdit.setHtml(QtWidgets.QApplication.translate("MainWindow", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'Segoe UI\'; font-size:10pt; font-weight:400; font-style:normal;\">\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'MS Shell Dlg 2\'; font-size:8.25pt;\">500</span></p></body></html>", None, -1))
        self.label_19.setText(QtWidgets.QApplication.translate("MainWindow", "Number of Steps: ", None, -1))
        self.groupBox_9.setTitle(QtWidgets.QApplication.translate("MainWindow", "Reporters", None, -1))
        self.State_Data_Reporter_Options_groupBox.setTitle(QtWidgets.QApplication.translate("MainWindow", "StateDataReporter Options", None, -1))
        self.label_26.setText(QtWidgets.QApplication.translate("MainWindow", "StateData Frequency: ", None, -1))
        self.StateData_frequency_textEdit.setHtml(QtWidgets.QApplication.translate("MainWindow", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'MS Shell Dlg 2\'; font-size:8.25pt; font-weight:400; font-style:normal;\">\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">100</p></body></html>", None, -1))
        self.DCD_Reporter_Options_groupBox.setTitle(QtWidgets.QApplication.translate("MainWindow", "DCD Reporter Options", None, -1))
        self.label_24.setText(QtWidgets.QApplication.translate("MainWindow", "DCD Write Frequency: ", None, -1))
        self.label_25.setText(QtWidgets.QApplication.translate("MainWindow", "DCD output: ", None, -1))
        self.DCD_write_freq_textEdit.setHtml(QtWidgets.QApplication.translate("MainWindow", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'Segoe UI\'; font-size:10pt; font-weight:400; font-style:normal;\">\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'MS Shell Dlg 2\'; font-size:8.25pt;\">100</span></p></body></html>", None, -1))
        self.DCD_Output_Name_textEdit.setHtml(QtWidgets.QApplication.translate("MainWindow", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'Segoe UI\'; font-size:10pt; font-weight:400; font-style:normal;\">\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'MS Shell Dlg 2\'; font-size:8.25pt;\">output.dcd</span></p></body></html>", None, -1))
        self.XTC_Reporter_Options_groupBox.setTitle(QtWidgets.QApplication.translate("MainWindow", "XTC Reporter Options", None, -1))
        self.label_36.setText(QtWidgets.QApplication.translate("MainWindow", "XTC Write Frequency: ", None, -1))
        self.label_37.setText(QtWidgets.QApplication.translate("MainWindow", "XTC output: ", None, -1))
        self.XTC_write_freq_textEdit.setHtml(QtWidgets.QApplication.translate("MainWindow", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'Segoe UI\'; font-size:10pt; font-weight:400; font-style:normal;\">\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'MS Shell Dlg 2\'; font-size:8.25pt;\">100</span></p></body></html>", None, -1))
        self.XTC_Output_Name_textEdit.setHtml(QtWidgets.QApplication.translate("MainWindow", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'Segoe UI\'; font-size:10pt; font-weight:400; font-style:normal;\">\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'MS Shell Dlg 2\'; font-size:8.25pt;\">output.xtc</span></p></body></html>", None, -1))
        self.label_22.setText(QtWidgets.QApplication.translate("MainWindow", "DCD Reporter: ", None, -1))
        self.label_40.setText(QtWidgets.QApplication.translate("MainWindow", "XTC Reporter: ", None, -1))
        self.label_23.setText(QtWidgets.QApplication.translate("MainWindow", "StateDataReporter: ", None, -1))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_4), QtWidgets.QApplication.translate("MainWindow", "Simulation", None, -1))
        self.Perturbation_TabWidget.setTabText(self.Perturbation_TabWidget.indexOf(self.Perturb_Advanced), QtWidgets.QApplication.translate("MainWindow", "Advanced", None, -1))
        self.visualization_settings_groupBox.setTitle(QtWidgets.QApplication.translate("MainWindow", "Visualization Settings", None, -1))
        self.activate_pymol_navigation.setText(QtWidgets.QApplication.translate("MainWindow", "Activate", None, -1))
        self.deactivate_pymol_navigation.setText(QtWidgets.QApplication.translate("MainWindow", "Deactivate", None, -1))
        self.refresh_pushButton.setText(QtWidgets.QApplication.translate("MainWindow", "Refresh", None, -1))
        self.ss_beatiful_snapshoot.setText(QtWidgets.QApplication.translate("MainWindow", "Beatiful Snap", None, -1))
        self.save_as_png_pushButton.setText(QtWidgets.QApplication.translate("MainWindow", "Save as png", None, -1))
        self.figure_settings_groupBox.setTitle(QtWidgets.QApplication.translate("MainWindow", "Figure Settings", None, -1))
        self.pymol_width_label.setText(QtWidgets.QApplication.translate("MainWindow", "Width: 1200", None, -1))
        self.pymol_height_label.setText(QtWidgets.QApplication.translate("MainWindow", "Height: 1080", None, -1))
        self.pymol_dpi_label.setText(QtWidgets.QApplication.translate("MainWindow", "Dpi: 150", None, -1))
        self.pymol_ray_label.setText(QtWidgets.QApplication.translate("MainWindow", "Ray: 1", None, -1))
        self.get_figure_pushButton.setText(QtWidgets.QApplication.translate("MainWindow", "Get Figure", None, -1))
        self.groupBox_2.setTitle(QtWidgets.QApplication.translate("MainWindow", "Each Residue Energy Calculation", None, -1))
        self.groupBox_29.setTitle(QtWidgets.QApplication.translate("MainWindow", "Start/ Stop", None, -1))
        self.stop_pushButton_3.setText(QtWidgets.QApplication.translate("MainWindow", "Stop", None, -1))
        self.network_calculate_pushButton.setText(QtWidgets.QApplication.translate("MainWindow", "CALCULATE", None, -1))
        self.label_100.setText(QtWidgets.QApplication.translate("MainWindow", "Threading:", None, -1))
        self.checkBox_5.setText(QtWidgets.QApplication.translate("MainWindow", "All CPU", None, -1))
        self.groupBox_31.setTitle(QtWidgets.QApplication.translate("MainWindow", "Residue Conservation", None, -1))
        self.label_102.setText(QtWidgets.QApplication.translate("MainWindow", "Threshold:", None, -1))
        self.label_133.setText(QtWidgets.QApplication.translate("MainWindow", "Chain ID:", None, -1))
        self.get_conserv_score_pushButton.setText(QtWidgets.QApplication.translate("MainWindow", "GET SCORES", None, -1))
        self.label_103.setText(QtWidgets.QApplication.translate("MainWindow", "PDB ID:", None, -1))
        self.label_104.setText(QtWidgets.QApplication.translate("MainWindow", "Conservation Scores", None, -1))
        self.residues_conservation_tableWidget.horizontalHeaderItem(0).setText(QtWidgets.QApplication.translate("MainWindow", "Residues", None, -1))
        self.residues_conservation_tableWidget.horizontalHeaderItem(1).setText(QtWidgets.QApplication.translate("MainWindow", "Scores", None, -1))
        self.use_conservation_checkBox.setText(QtWidgets.QApplication.translate("MainWindow", "Use Conservation", None, -1))
        self.groupBox_30.setTitle(QtWidgets.QApplication.translate("MainWindow", "Network Calculation Parameters", None, -1))
        self.label_98.setText(QtWidgets.QApplication.translate("MainWindow", "Node Threshold:", None, -1))
        self.label_99.setText(QtWidgets.QApplication.translate("MainWindow", "CutOff", None, -1))
        self.node_threshold_checkBox.setText(QtWidgets.QApplication.translate("MainWindow", "None", None, -1))
        self.label_105.setText(QtWidgets.QApplication.translate("MainWindow", "Source Residue", None, -1))
        self.label_97.setText(QtWidgets.QApplication.translate("MainWindow", "Target Residu(s)", None, -1))
        self.label_101.setText(QtWidgets.QApplication.translate("MainWindow", "Selected Residue(s)", None, -1))
        self.groupBox_28.setTitle(QtWidgets.QApplication.translate("MainWindow", "Output", None, -1))
        self.label_94.setText(QtWidgets.QApplication.translate("MainWindow", "Output Directory:", None, -1))
        self.label_95.setText(QtWidgets.QApplication.translate("MainWindow", "PPI Network Name:", None, -1))
        self.output_directory_button.setText(QtWidgets.QApplication.translate("MainWindow", "Browse", None, -1))
        self.PPI_Network_name_lineedit.setText(QtWidgets.QApplication.translate("MainWindow", "protein_general_cutoff_network.gml", None, -1))
        self.groupBox_27.setTitle(QtWidgets.QApplication.translate("MainWindow", "Inputs", None, -1))
        self.label_93.setText(QtWidgets.QApplication.translate("MainWindow", "Response Time File:", None, -1))
        self.label_91.setText(QtWidgets.QApplication.translate("MainWindow", "Protein Bound Form:", None, -1))
        self.upload_boundForm_pdb_Button.setText(QtWidgets.QApplication.translate("MainWindow", "Upload", None, -1))
        self.response_time_upload_Button.setText(QtWidgets.QApplication.translate("MainWindow", "Upload", None, -1))
        self.add_tabb.setText(QtWidgets.QApplication.translate("MainWindow", "PushButton", None, -1))
        self.analysis_TabWidget.setTabText(self.analysis_TabWidget.indexOf(self.Perturb_Quick_3), QtWidgets.QApplication.translate("MainWindow", "Inputs", None, -1))
        self.label_28.setText(QtWidgets.QApplication.translate("MainWindow", "Intersection Shortest Path(s)", None, -1))
        self.label_29.setText(QtWidgets.QApplication.translate("MainWindow", "Shortest Path(s)", None, -1))
        self.analysis_TabWidget.setTabText(self.analysis_TabWidget.indexOf(self.tab_5), QtWidgets.QApplication.translate("MainWindow", "Page", None, -1))
        self.label_34.setText(QtWidgets.QApplication.translate("MainWindow", "ABOUT & CONTACT", None, -1))
        self.label_33.setText(QtWidgets.QApplication.translate("MainWindow", "MARMARA UNIVERSITY", None, -1))
        self.label_41.setText(QtWidgets.QApplication.translate("MainWindow", "COMPUTATIONAL BIOLOGY AND BIOINFORMATICS RESEARCH GROUP", None, -1))
        self.textBrowser.setHtml(QtWidgets.QApplication.translate("MainWindow", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'Segoe UI\'; font-size:8.25pt; font-weight:400; font-style:normal;\">\n"
"<table border=\"0\" style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px;\" cellspacing=\"2\" cellpadding=\"0\">\n"
"<tr>\n"
"<td>\n"
"<p align=\"justify\" style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-family:\'MS Shell Dlg 2\'; font-size:8pt; font-weight:600; color:#ffffff;\"><br /></p>\n"
"<p align=\"justify\" style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><a name=\"LC1\"></a><span style=\" font-family:\'MS Shell Dlg 2\'; font-size:8pt; font-weight:600; color:#ffffff;\">M</span><span style=\" font-family:\'MS Shell Dlg 2\'; font-size:8pt; font-weight:600; color:#ffffff;\">IT License</span></p></td></tr>\n"
"<tr>\n"
"<td></td></tr></table>\n"
"<table border=\"0\" style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px;\" cellspacing=\"2\" cellpadding=\"0\">\n"
"<tr>\n"
"<td>\n"
"<p align=\"justify\" style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><a name=\"LC2\"></a><span style=\" font-family:\'MS Shell Dlg 2\'; font-size:8pt; color:#ffffff;\">C</span><span style=\" font-family:\'MS Shell Dlg 2\'; font-size:8pt; color:#ffffff;\">opyright (c) 2020 </span><span style=\" font-family:\'MS Shell Dlg 2\'; font-size:8pt; font-weight:600; color:#ffffff;\">Halil İbrahim ÖZDEMİR</span></p></td></tr>\n"
"<tr>\n"
"<td></td></tr></table>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-size:8pt;\"><br /></p>\n"
"<table border=\"0\" style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px;\" cellspacing=\"2\" cellpadding=\"0\">\n"
"<tr>\n"
"<td>\n"
"<p align=\"justify\" style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><a name=\"LC4\"></a><span style=\" font-family:\'MS Shell Dlg 2\'; font-size:8pt; color:#ffffff;\">P</span><span style=\" font-family:\'MS Shell Dlg 2\'; font-size:8pt; color:#ffffff;\">ermission is hereby granted, free of charge, to any person obtaining a copy </span><a name=\"LC6\"></a><span style=\" font-family:\'MS Shell Dlg 2\'; font-size:8pt; color:#ffffff;\">o</span><span style=\" font-family:\'MS Shell Dlg 2\'; font-size:8pt; color:#ffffff;\">f this software and associated documentation files (the &quot;Software&quot;), to deal </span><a name=\"LC7\"></a><span style=\" font-family:\'MS Shell Dlg 2\'; font-size:8pt; color:#ffffff;\">i</span><span style=\" font-family:\'MS Shell Dlg 2\'; font-size:8pt; color:#ffffff;\">n the Software without restriction, including without limitation the rights </span><a name=\"LC8\"></a><span style=\" font-family:\'MS Shell Dlg 2\'; font-size:8pt; color:#ffffff;\">t</span><span style=\" font-family:\'MS Shell Dlg 2\'; font-size:8pt; color:#ffffff;\">o use, copy, modify, merge, publish, distribute, sublicense, and/or sell </span><a name=\"LC9\"></a><span style=\" font-family:\'MS Shell Dlg 2\'; font-size:8pt; color:#ffffff;\">c</span><span style=\" font-family:\'MS Shell Dlg 2\'; font-size:8pt; color:#ffffff;\">opies of the Software, and to permit persons to whom the Software is </span><a name=\"LC10\"></a><span style=\" font-family:\'MS Shell Dlg 2\'; font-size:8pt; color:#ffffff;\">f</span><span style=\" font-family:\'MS Shell Dlg 2\'; font-size:8pt; color:#ffffff;\">urnished to do so, subject to the following conditions: </span><a name=\"LC11\"></a><span style=\" font-family:\'MS Shell Dlg 2\'; font-size:8pt; color:#ffffff;\">T</span><span style=\" font-family:\'MS Shell Dlg 2\'; font-size:8pt; color:#ffffff;\">he above copyright notice and this permission notice shall be included in all </span><a name=\"LC13\"></a><span style=\" font-family:\'MS Shell Dlg 2\'; font-size:8pt; color:#ffffff;\">c</span><span style=\" font-family:\'MS Shell Dlg 2\'; font-size:8pt; color:#ffffff;\">opies or substantial portions of the Software.</span></p></td></tr>\n"
"<tr>\n"
"<td></td></tr></table>\n"
"<table border=\"0\" style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px;\" cellspacing=\"2\" cellpadding=\"0\">\n"
"<tr>\n"
"<td>\n"
"<p align=\"justify\" style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><a name=\"LC14\"></a><span style=\" font-family:\'MS Shell Dlg 2\'; font-size:8pt; font-weight:600; color:#ffffff;\">T</span><span style=\" font-family:\'MS Shell Dlg 2\'; font-size:8pt; font-weight:600; color:#ffffff;\">HE SOFTWARE IS PROVIDED &quot;AS IS&quot;, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR </span><a name=\"LC16\"></a><span style=\" font-family:\'MS Shell Dlg 2\'; font-size:8pt; font-weight:600; color:#ffffff;\">I</span><span style=\" font-family:\'MS Shell Dlg 2\'; font-size:8pt; font-weight:600; color:#ffffff;\">MPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, </span><a name=\"LC17\"></a><span style=\" font-family:\'MS Shell Dlg 2\'; font-size:8pt; font-weight:600; color:#ffffff;\">F</span><span style=\" font-family:\'MS Shell Dlg 2\'; font-size:8pt; font-weight:600; color:#ffffff;\">ITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE </span><a name=\"LC18\"></a><span style=\" font-family:\'MS Shell Dlg 2\'; font-size:8pt; font-weight:600; color:#ffffff;\">A</span><span style=\" font-family:\'MS Shell Dlg 2\'; font-size:8pt; font-weight:600; color:#ffffff;\">UTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER </span><a name=\"LC19\"></a><span style=\" font-family:\'MS Shell Dlg 2\'; font-size:8pt; font-weight:600; color:#ffffff;\">L</span><span style=\" font-family:\'MS Shell Dlg 2\'; font-size:8pt; font-weight:600; color:#ffffff;\">IABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, </span><a name=\"LC19\"></a><span style=\" font-family:\'MS Shell Dlg 2\'; font-size:8pt; font-weight:600; color:#ffffff;\">L</span><span style=\" font-family:\'MS Shell Dlg 2\'; font-size:8pt; font-weight:600; color:#ffffff;\">IABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, </span><a name=\"LC20\"></a><span style=\" font-family:\'MS Shell Dlg 2\'; font-size:8pt; font-weight:600; color:#ffffff;\">O</span><span style=\" font-family:\'MS Shell Dlg 2\'; font-size:8pt; font-weight:600; color:#ffffff;\">UT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.</span></p>\n"
"<table border=\"0\" style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px;\" cellspacing=\"2\" cellpadding=\"0\">\n"
"<tr>\n"
"<td>\n"
"<table border=\"0\" style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px;\" cellspacing=\"2\" cellpadding=\"0\">\n"
"<tr>\n"
"<td></td></tr>\n"
"<tr>\n"
"<td></td></tr></table></td></tr>\n"
"<tr>\n"
"<td></td></tr></table></td></tr>\n"
"<tr>\n"
"<td></td></tr></table></body></html>", None, -1))
        self.commandLinkButton_3.setText(QtWidgets.QApplication.translate("MainWindow", "GO TO PERSONAL WEB PAGE", None, -1))
        self.commandLinkButton_2.setText(QtWidgets.QApplication.translate("MainWindow", "GO TO PERSONAL WEB PAGE", None, -1))
        self.commandLinkButton.setText(QtWidgets.QApplication.translate("MainWindow", "CREDITS", None, -1))
        self.label_credits.setText(QtWidgets.QApplication.translate("MainWindow", "Marmara University / Bioengineering / Computational Biology and Bioinformatics", None, -1))
        self.label_version.setText(QtWidgets.QApplication.translate("MainWindow", "Version 0.1", None, -1))

import icons_path_rc
import files_rc
