## ==> GUI FILE
from gui.ui_styles import Style

## ==> GLOBALS
GLOBAL_STATE = 0
GLOBAL_TITLE_BAR = True

## ==> COUT INITIAL MENU
count = 1
from .app_functions import *
from .PyMolWidget import PymolQtWidget
from analysis import VisJS_Widget
from .message import Message_Boxes


class UIFunctions(MainWindow):
    # ----- > GLOBALS
    GLOBAL_STATE = 0
    GLOBAL_TITLE_BAR = True
    ProteinView = False
    Protein3DNetworkView = False
    VisJSEngineView = False

    ####################################################################################################################
    #                                       == > START - GUI FUNCTIONS < ==                                           #
    ####################################################################################################################

    # ----- > Close Application
    @staticmethod
    def close_application(self):
        try:
            """Close Application Question Message Box."""
            close_answer = Message_Boxes.Question_message(self, "Be carefull!", "Do you really want to close the "
                                                                                "program?", Style.MessageBox_stylesheet)

            if close_answer == QMessageBox.Yes:
                self.thread_main.exit()
                self.close()
            if close_answer == QMessageBox.No:
                pass

        except Exception as inst:
            Message_Boxes.Critical_message(self, 'An unexpected error has occurred!', str(inst),
                                           Style.MessageBox_stylesheet)

    # ----- > Maximize / Restore
    def maximize_restore(self):
        global GLOBAL_STATE
        status = GLOBAL_STATE

        if self.isMaximized():
            GLOBAL_STATE = 0
            self.showNormal()
            # self.resize(self.width() + 1, self.height() + 1)
            # self.horizontalLayout.setContentsMargins(10, 10, 10, 10)
            self.btn_maximize_restore.setToolTip("Maximize")
            self.btn_maximize_restore.setIcon(QtGui.QIcon(u":/16x16/icons/16x16/cil-window-maximize.png"))
            # self.frame_top_btns.setStyleSheet("background-color: rgba(27, 29, 35, 200)")
            self.frame_size_grip.show()

        else:
            self.showMaximized()
            GLOBAL_STATE = 1
            self.horizontalLayout.setContentsMargins(0, 0, 0, 0)
            self.btn_maximize_restore.setToolTip("Restore")
            self.btn_maximize_restore.setIcon(QtGui.QIcon(u":/16x16/icons/16x16/cil-window-restore.png"))
            # self.frame_top_btns.setStyleSheet("background-color: rgb(27, 29, 35)")
            self.frame_size_grip.hide()

    # ----- > Return Status
    def returStatus():
        return GLOBAL_STATE

    # ----- > Set Status
    def setStatus(status):
        global GLOBAL_STATE
        GLOBAL_STATE = status

    # ----- > Enable Maximum Size
    def enableMaximumSize(self, width, height):
        if width != '' and height != '':
            self.setMaximumSize(QSize(width, height))
            self.frame_size_grip.hide()
            self.btn_maximize_restore.hide()

    # ----- > Toggle Menu
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

    # ----- > Set Title Bar
    def removeTitleBar(self):
        global GLOBAL_TITLE_BAR
        GLOBAL_TITLE_BAR = self

    # ----- > Header Texts
    # ==> LABEL TITLE
    def labelTitle(self, text):
        self.label_title_bar_top.setText(text)

    # ==> LABEL DESCRIPTION
    def labelDescription(self, text):
        self.label_top_info_1.setText(text)

    # ----- > Dynamic Menus
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

    # ----- > Select Menu
    def selectMenu(getStyle):
        select = getStyle + ("QPushButton { border-right: 7px solid rgb(44, 49, 60); }")
        return select

    # ----- > Deselect Menu
    def deselectMenu(getStyle):
        deselect = getStyle.replace("QPushButton { border-right: 7px solid rgb(44, 49, 60); }", "")
        return deselect

    # ----- > Start Selection
    def selectStandardMenu(self, widget):
        for w in self.frame_left_menu.findChildren(QPushButton):
            if w.objectName() == widget:
                w.setStyleSheet(UIFunctions.selectMenu(w.styleSheet()))

    # ----- > Reset Selection
    def resetStyle(self, widget):
        for w in self.frame_left_menu.findChildren(QPushButton):
            if w.objectName() != widget:
                w.setStyleSheet(UIFunctions.deselectMenu(w.styleSheet()))

    # ----- > Change Page Label Text
    def labelPage(self, text):
        newText = '| ' + text.upper()
        self.label_top_info_2.setText(newText)

    # ----- > User Icon
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

    ####################################################################################################################
    #                                       == > END - GUI FUNCTIONS < ==                                           #
    ####################################################################################################################
    ####################################################################################################################
    #                                       == > START - GUI DEFINITIONS < ==                                          #
    ####################################################################################################################

    # ----- > UI Definitions
    def uiDefinitions(self):
        def dobleClickMaximizeRestore(event):
            # IF DOUBLE CLICK CHANGE STATUS
            if event.type() == QtCore.QEvent.MouseButtonDblClick:
                QtCore.QTimer.singleShot(250, lambda: UIFunctions.maximize_restore(self))

        # REMOVE ==> STANDARD TITLE BAR
        if GLOBAL_TITLE_BAR:
            self.setStyleSheet("background:rgb(27, 29, 35);")
            self.setWindowFlags(QtCore.Qt.FramelessWindowHint)
            self.frame_main.setAttribute(QtCore.Qt.WA_TranslucentBackground)
            self.centralwidget.setAttribute(QtCore.Qt.WA_NoSystemBackground)
            # self.setStyleSheet("background:rgb(27, 29, 35);")
            self.frame_label_top_btns.mouseDoubleClickEvent = dobleClickMaximizeRestore
        else:
            self.horizontalLayout.setContentsMargins(0, 0, 0, 0)
            self.frame_label_top_btns.setContentsMargins(8, 0, 0, 5)
            self.frame_label_top_btns.setMinimumHeight(42)
            self.frame_icon_top_bar.hide()
            self.frame_btns_right.hide()
            self.frame_size_grip.hide()

        # SHOW ==> DROP SHADOW
        # self.shadow = QGraphicsDropShadowEffect(self)
        # self.shadow.setBlurRadius(0)
        # self.shadow.setXOffset(0)
        # self.shadow.setYOffset(0)
        # self.shadow.setColor(QColor(0, 0, 0, 50))
        # self.frame_main.setGraphicsEffect(self.shadow)
        # DROP SHADOW EFFECT

        ## ==> RESIZE WINDOW
        self.sizegrip = QSizeGrip(self.frame_size_grip)
        self.sizegrip.setStyleSheet("width: 20px; height: 20px; margin 0px; padding: 0px;")

        ### ==> MINIMIZE
        self.btn_minimize.clicked.connect(lambda: self.showMinimized())

        ## ==> MAXIMIZE/RESTORE
        self.btn_maximize_restore.clicked.connect(lambda: UIFunctions.maximize_restore(self))

        ## SHOW ==> CLOSE APPLICATION
        self.btn_close.clicked.connect(lambda: UIFunctions.close_application(self))

    ####################################################################################################################
    #                                         == > END - GUI DEFINITIONS < ==                                          #
    ####################################################################################################################
    ####################################################################################################################
    #                             == > START - OPEN SOURCE PYMOL 2.4 INTEGRATION < ==                                  #
    ####################################################################################################################
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
            self.ProteinView.initial_pymol_visual()

            self.Protein3DNetworkView = PymolQtWidget(self)
            verticalLayoutProteinNetworkView = QVBoxLayout(self.PyMOL_3D_network_Widget)
            verticalLayoutProteinNetworkView.addWidget(self.Protein3DNetworkView)
            self.setLayout(verticalLayoutProteinNetworkView)
            self.Protein3DNetworkView.update()
            self.Protein3DNetworkView.show()
            verticalLayoutProteinNetworkView.setContentsMargins(0, 0, 0, 0)
            # self.Protein3DNetworkView.initial_pymol_visual()

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

    def load_pdb_to_3DNetwork(self, pdb_file):
        self.Protein3DNetworkView.reinitialize()
        self.Protein3DNetworkView.loadMolFile(pdb_file)
        self.Protein3DNetworkView.update()
        self.Protein3DNetworkView.show()

    def start_VisJS_2D_Network(self):
        self.VisJSEngineView = VisJS_Widget.VisJS_QtWidget()
        self.Network_2D_verticalLayout.addWidget(self.VisJSEngineView.m_output)

    # def load_nx_to_VisJS_2D_Network(self, intersection_graph='2d_network.html', gml_file='example.gml'):
    #     initial_2d_network_html_directory = os.path.join(os.getcwd(), 'analysis')
    #     initial_2d_network_html_path = os.path.join(initial_2d_network_html_directory, intersection_graph)
    #
    #     gml_graph = nx.read_gml(os.path.join(initial_2d_network_html_directory, gml_file))
    #     self.VisJSEngineView.load_network_component(network=gml_graph, html_file=initial_2d_network_html_path)

    def load_pdb_to_pymol(self, pdb_file):
        self.ProteinView.reinitialize()
        self.ProteinView.loadMolFile(pdb_file)
        self.ProteinView.update()
        self.ProteinView.show()

    def show_residue_labels(self):
        itemsTextList = self.selected_residues_listWidget.currentItem().text()
        # [str(self.selected_residues_listWidget.item(i).text()) for i in
        #              range(self.selected_residues_listWidget.count())]

        self.ProteinView.selection_color(itemsTextList)
        self.ProteinView.update()

    def clear_residue_labels(self):
        self.ProteinView.clear_all_labels()
        self.ProteinView.update()

    def activate_navigation_on_Pymol(self):
        self.ProteinView.activate_navigation_tool()
        self.ProteinView.paintGL()
        self.ProteinView.update()
        self.ProteinView.show()

    def deactivate_navigation_on_Pymol(self):
        self.ProteinView.deactivate_navigation_tool()
        self.ProteinView.paintGL()
        self.ProteinView.update()
        self.ProteinView.show()

    def show_beatiful_in_Pymol(self):
        self.ProteinView.set_ss_figure()
        self.ProteinView.update()
        self.ProteinView.show()

    def save_as_png_Pymol(self):
        filedialog = QFileDialog(self)
        filedialog.setDefaultSuffix("png")
        filedialog.setNameFilter("PNG Files (*.png);;All files (*.*)")
        filedialog.setAcceptMode(QFileDialog.AcceptSave)
        selected = filedialog.exec()

        if selected:
            filename = filedialog.selectedFiles()[0]
        else:
            return
        if filename == "":
            Message_Boxes.Warning_message(self, 'png save failed!', "No file name selected.",
                                          Style.MessageBox_stylesheet)
            return

        try:
            self.ProteinView.get_png_figure(filename, width=self.width_horizontalSlider.value(),
                                            height=self.height_horizontalSlider.value(),
                                            dpi=self.dpi_horizontalSlider.value(),
                                            ray=self.ray_horizontalSlider.value())
            self.ProteinView.update()

        except Exception as save_err:
            Message_Boxes.Critical_message(self, 'png save failed!', str(save_err), Style.MessageBox_stylesheet)

    """
    def deleteLayout(self, verticalLayoutProteinView):
        if verticalLayoutProteinView is not None:
            sip.delete(verticalLayoutProteinView)
    """
    ####################################################################################################################
    #                               == > END - OPEN SOURCE PYMOL 2.4 INTEGRATION < ==                                  #
    ####################################################################################################################
