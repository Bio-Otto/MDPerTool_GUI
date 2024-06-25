# class NetworkTools:
#     def plot_networks(self):
#         clean_graph_list = self._filter_graphs_by_node_threshold()
#         clean_log_list = self._filter_logs_by_node_threshold()
#
#         if clean_graph_list:
#             self._create_intersection_graph(clean_graph_list)
#             self._setup_analysis_tab(clean_graph_list, clean_log_list)
#
#     def _filter_graphs_by_node_threshold(self):
#         clean_graph_list = []
#         if self.node_threshold is not None:
#             for graph in self.network_holder:
#                 if len(graph.nodes()) > self.node_threshold:
#                     clean_graph_list.append(graph)
#         else:
#             for graph in self.network_holder:
#                 if len(graph.nodes()) > 0:
#                     clean_graph_list.append(graph)
#         return clean_graph_list
#
#     def _filter_logs_by_node_threshold(self):
#         clean_log_list = []
#         if self.node_threshold is not None:
#             for log, graph in zip(self.log_holder, self.network_holder):
#                 if len(graph.nodes()) > self.node_threshold:
#                     clean_log_list.append(log)
#         else:
#             for log, graph in zip(self.log_holder, self.network_holder):
#                 if len(graph.nodes()) > 0:
#                     clean_log_list.append(log)
#         return clean_log_list
#
#     def _create_intersection_graph(self, graph_list):
#         # Implementation for creating intersection graph and writing to GML file goes here
#         pass
#
#     def _setup_analysis_tab(self, graph_list, log_list):
#         import PyQt5.QtWidgets as QtWidgets
#         import PyQt5.QtCore as QtCore
#
#         tab = QtWidgets.QWidget()
#         tab.setObjectName(f"Analysis_{self.tab_count_on_analysis}")
#
#         self.analysis_TabWidget.tabBar().setTabButton(0, QtWidgets.QTabBar.RightSide, None)
#         self.tab_count_on_analysis = self.analysis_TabWidget.count()
#
#         horizontalLayout = QtWidgets.QHBoxLayout(tab)
#         horizontalLayout.setObjectName(f"horizontalLayout_{self.tab_count_on_analysis}")
#         gridLayout = QtWidgets.QGridLayout()
#         gridLayout.setObjectName(f"gridLayout_{self.tab_count_on_analysis}")
#
#         self._add_shortest_path_list(gridLayout)
#         self._add_labels(gridLayout)
#         self._add_dissipation_widget(gridLayout)
#         self._add_navigation_buttons(gridLayout, tab)
#         self._add_analysis_settings_groupbox(gridLayout, tab)
#
#         horizontalLayout.addLayout(gridLayout)
#         self.analysis_TabWidget.addTab(tab, f"Analysis {self.tab_count_on_analysis}")
#
#     def _add_shortest_path_list(self, layout):
#
#         shortest_path_listWidget = QtWidgets.QListWidget()
#         shortest_path_listWidget.setMaximumSize(QtCore.QSize(450, 16777215))
#         shortest_path_listWidget.setObjectName("shortest_path_listWidget")
#         layout.addWidget(shortest_path_listWidget, 1, 0, 1, 1)
#
#     def _add_labels(self, layout):
#         import PyQt5.QtWidgets as QtWidgets
#         import PyQt5.QtCore as QtCore
#
#         label_stylesheet = (
#             "QLabel {"
#             "    background-color: rgb(27, 29, 35);"
#             "    border-radius: 5px;"
#             "    border: 2px solid rgb(27, 29, 35);"
#             "    padding: 1px 1px 1px 1px;"
#             "    border-bottom-color: rgb(157, 90, 198);"
#             "}"
#             "QLabel:hover{"
#             "    border: 2px solid rgb(64, 71, 88);"
#             "    selection-color: rgb(127, 5, 64);"
#             "}"
#         )
#
#         label_shortest_path = QtWidgets.QLabel("Shortest Path(s)")
#         label_shortest_path.setMinimumSize(QtCore.QSize(0, 22))
#         label_shortest_path.setMaximumSize(QtCore.QSize(450, 22))
#         label_shortest_path.setStyleSheet(label_stylesheet)
#         label_shortest_path.setObjectName("label_shortest_path")
#         layout.addWidget(label_shortest_path, 0, 0, 1, 1)
#
#         label_energy_dissipation = QtWidgets.QLabel("Energy Dissipation Curve")
#         label_energy_dissipation.setMinimumSize(QtCore.QSize(0, 22))
#         label_energy_dissipation.setMaximumSize(QtCore.QSize(450, 22))
#         label_energy_dissipation.setStyleSheet(label_stylesheet)
#         label_energy_dissipation.setObjectName("label_energy_dissipation")
#         layout.addWidget(label_energy_dissipation, 4, 0, 1, 1)
#
#     def _add_dissipation_widget(self, layout):
#
#         dissipation_curve_widget = WidgetPlot(self)
#
#         toolbarSizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
#         dissipation_curve_widget.toolbar.setSizePolicy(toolbarSizePolicy)
#
#         dissipationVerticalLayout = QtWidgets.QVBoxLayout()
#         dissipationVerticalLayout.addWidget(dissipation_curve_widget.toolbar)
#         dissipationVerticalLayout.addWidget(dissipation_curve_widget.canvas)
#         dissipation_curve_widget.setLayout(dissipationVerticalLayout)
#
#         sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
#         dissipation_curve_widget.setSizePolicy(sizePolicy)
#         dissipation_curve_widget.setMinimumSize(QtCore.QSize(0, 400))
#         dissipation_curve_widget.setMaximumSize(QtCore.QSize(450, 450))
#         dissipation_curve_widget.setObjectName("dissipation_curve_widget")
#         layout.addWidget(dissipation_curve_widget, 5, 0, 1, 2)
#
#     def _add_navigation_buttons(self, layout, tab):
#
#         button_stylesheet = (
#             "QPushButton {"
#             "    color: white;"
#             "    border: 2px solid rgb(52, 59, 72);"
#             "    border-radius: 5px;"
#             "    background-color: rgb(110, 105, 225);"
#             "    border-width: 1px;"
#             "    outline: none;"
#             "}"
#             "QPushButton:hover {"
#             "    background-color: rgb(22, 200, 244);"
#             "    border: 2px solid rgb(61, 70, 86);"
#             "}"
#             "QPushButton:pressed {"
#             "    background-color: rgb(15, 133, 163);"
#             "    border: 2px solid rgb(43, 50, 61);"
#             "}"
#         )
#
#         show_navigation_button = QtWidgets.QPushButton(tab)
#         show_navigation_button.setSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Preferred)
#         show_navigation_button.setMaximumSize(QtCore.QSize(20, 61))
#         show_navigation_button.setStyleSheet(button_stylesheet)
#         show_navigation_button.setIcon(QtGui.QIcon(":/16x16/icons/16x16/cil-chevron-circle-right-alt.png"))
#         show_navigation_button.setObjectName("show_navigation_button")
#         layout.addWidget(show_navigation_button, 0, 1, 6, 1)
#
#         hide_navigation_button = QtWidgets.QPushButton(tab)
#         hide_navigation_button.setMaximumSize(QtCore.QSize(20, 61))
#         hide_navigation_button.setStyleSheet(button_stylesheet)
#         hide_navigation_button.setIcon(QtGui.QIcon(":/16x16/icons/16x16/cil-chevron-circle-right-alt.png"))
#         hide_navigation_button.setObjectName("hide_navigation")
#         layout.addWidget(hide_navigation_button, 0, 3, 6, 1)
#
#
# def _add_analysis_settings_groupbox(self, layout, tab):
#
#     groupbox_stylesheet = (
#         "QGroupBox {"
#         "    border: 1px solid black;"
#         "    border-radius: 5px;"
#         "    border-top-color: rgb(157, 90, 198);"
#         "    border-left-color: rgb(157, 90, 198);"
#         "    border-bottom-color: rgb(157, 90, 198);"
#         "    border-right-color: rgb(157, 90, 198);"
#         "}"
#     )
#
#     button_stylesheet = (
#         "QPushButton {"
#         "    color: white;"
#         "    border: 2px solid rgb(52, 59, 72);"
#         "    border-radius: 5px;"
#         "    background-color: rgb(22, 200, 244);"
#         "    margin-top: 1px;"
#         "    margin-bottom: 1px;"
#         "    border-width: 1px;"
#         "    padding: 5px;"
#         "    outline: none;"
#         "}"
#         "QPushButton:hover {"
#         "    background-color: rgb(255, 17, 100);"
#         "    border: 2px solid rgb(61, 70, 86);"
#         "}"
#         "QPushButton:pressed {"
#         "    background-color: rgb(15, 133, 163);"
#         "    border: 2px solid rgb(43, 50, 61);"
#         "}"
#     )
#
#     analysis_settings_groupBox = QtWidgets.QGroupBox(tab)
#     analysis_settings_groupBox.setTitle('Visualization Settings')
#     analysis_settings_groupBox.setMinimumSize(QtCore.QSize(170, 0))
#     analysis_settings_groupBox.setMaximumSize(QtCore.QSize(170, 16777215))
#     analysis_settings_groupBox.setStyleSheet(groupbox_stylesheet)
#     analysis_settings_groupBox.setAlignment(QtCore.Qt.AlignLeading | QtCore.Qt.AlignLeft | QtCore.Qt.AlignVCenter)
#     analysis_settings_groupBox.setObjectName("analysis_settings_groupBox")
#
#     gridLayoutWidget_on_analysis = QtWidgets.QWidget(analysis_settings_groupBox)
#     gridLayoutWidget_on_analysis.setGeometry(QtCore.QRect(11, 50, 151, 251))
#     gridLayoutWidget_on_analysis.setObjectName("gridLayoutWidget_on_analysis")
#     verticalLayout_analysis = QtWidgets.QVBoxLayout(gridLayoutWidget_on_analysis)
#     verticalLayout_analysis.setContentsMargins(0, 0, 0, 0)
#     verticalLayout_analysis.setObjectName("verticalLayout_analysis")
#
#     activate_pymol_navigation_on_analysis = QtWidgets.QPushButton(gridLayoutWidget_on_analysis)
#     activate_pymol_navigation_on_analysis.setText('Activate')
#     activate_pymol_navigation_on_analysis.setSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
#     activate_pymol_navigation_on_analysis.setMinimumSize(QtCore.QSize(149, 33))
#     activate_pymol_navigation_on_analysis.setMaximumSize(QtCore.QSize(110, 33))
#     activate_pymol_navigation_on_analysis.setFont(QtGui.QFont("Segoe UI", 9))
#     activate_pymol_navigation_on_analysis.setStyleSheet(button_stylesheet)
#     activate_pymol_navigation_on_analysis.setIcon(QtGui.QIcon(":/24x24/icons/24x24/cil-cursor.png"))
#     activate_pymol_navigation_on_analysis.setObjectName("activate_pymol_navigation_on_analysis")
#     verticalLayout_analysis.addWidget(activate_pymol_navigation_on_analysis)
#
#     deactivate_pymol_navigation_on_analysis = QtWidgets.QPushButton(gridLayoutWidget_on_analysis)
#     deactivate_pymol_navigation_on_analysis.setText('Deactivate')
#     deactivate_pymol_navigation_on_analysis.setSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
#     deactivate_pymol_navigation_on_analysis.setMinimumSize(QtCore.QSize(149, 33))
#     deactivate_pymol_navigation_on_analysis.setMaximumSize(QtCore.QSize(110, 33))
#     deactivate_pymol_navigation_on_analysis.setFont(QtGui.QFont("Segoe UI", 9))
#     deactivate_pymol_navigation_on_analysis.setStyleSheet(button_stylesheet)
#     deactivate_pymol_navigation_on_analysis.setIcon(QtGui.QIcon(":/20x20/icons/20x20/cil-x.png"))
#     deactivate_pymol_navigation_on_analysis.setObjectName("deactivate_pymol_navigation")
#     verticalLayout_analysis.addWidget(deactivate_pymol_navigation_on_analysis)
#
#     refresh_pushButton_on_analysis = QtWidgets.QPushButton(gridLayoutWidget_on_analysis)
#     refresh_pushButton_on_analysis.setText('Refresh')
#     refresh_pushButton_on_analysis.setSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
#     refresh_pushButton_on_analysis.setMinimumSize(QtCore.QSize(149, 33))
#     refresh_pushButton_on_analysis.setMaximumSize(QtCore.QSize(110, 33))
#     refresh_pushButton_on_analysis.setFont(QtGui.QFont("Segoe UI", 9))
#     refresh_pushButton_on_analysis.setStyleSheet(button_stylesheet)
#     refresh_pushButton_on_analysis.setIcon(QtGui.QIcon(":/16x16/icons/16x16/cil-reload.png"))
#     refresh_pushButton_on_analysis.setObjectName("refresh_pushButton_on_analysis")
#     verticalLayout_analysis.addWidget(refresh_pushButton_on_analysis)
#
#     ss_beatiful_snapshoot_on_analysis = QtWidgets.QPushButton(gridLayoutWidget_on_analysis)
#     ss_beatiful_snapshoot_on_analysis.setText('Beautiful Snap')
#     ss_beatiful_snapshoot_on_analysis.setSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
#     ss_beatiful_snapshoot_on_analysis.setMinimumSize(QtCore.QSize(149, 33))
#     ss_beatiful_snapshoot_on_analysis.setMaximumSize(QtCore.QSize(110, 33))
#     ss_beatiful_snapshoot_on_analysis.setFont(QtGui.QFont("Segoe UI", 9))
#     ss_beatiful_snapshoot_on_analysis.setStyleSheet(button_stylesheet)
#     ss_beatiful_snapshoot_on_analysis.setIcon(QtGui.QIcon(":/16x16/icons/16x16/cil-camera.png"))
#     ss_beatiful_snapshoot_on_analysis.setObjectName("ss_beatiful_snapshoot_on_analysis")
#     verticalLayout_analysis.addWidget(ss_beatiful_snapshoot_on_analysis)
#
#     save_as_png_on_analysis_pushButton = QtWidgets.QPushButton(gridLayoutWidget_on_analysis)
#     save_as_png_on_analysis_pushButton.setText('Save as png')
#     save_as_png_on_analysis_pushButton.setSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
#     save_as_png_on_analysis_pushButton.setMinimumSize(QtCore.QSize(149, 33))
#     save_as_png_on_analysis_pushButton.setMaximumSize(QtCore.QSize(110, 33))
#     save_as_png_on_analysis_pushButton.setFont(QtGui.QFont("Segoe UI", 9))
#     save_as_png_on_analysis_pushButton.setStyleSheet(button_stylesheet)
#     save_as_png_on_analysis_pushButton.setIcon(QtGui.QIcon(":/16x16/icons/16x16/cil-data-transfer-down.png"))
#     save_as_png_on_analysis_pushButton.setObjectName("save_as_png_on_analysis_pushButton")
#     verticalLayout_analysis.addWidget(save_as_png_on_analysis_pushButton)
#
#     figure_settings_on_analysis_groupBox = QtWidgets.QGroupBox(analysis_settings_groupBox)
#     figure_settings_on_analysis_groupBox.setTitle('Figure Settings')
#     figure_settings_on_analysis_groupBox.setGeometry(QtCore.QRect(-1, 300, 171, 401))
#     figure_settings_on_analysis_groupBox.setObjectName("figure_settings_on_analysis_groupBox")
#
#     verticalLayoutWidget_on_analysis_2 = QtWidgets.QWidget(figure_settings_on_analysis_groupBox)
#     verticalLayoutWidget_on_analysis_2.setGeometry(QtCore.QRect(10, 54, 151, 311))
#     verticalLayoutWidget_on_analysis_2.setObjectName("verticalLayoutWidget_on_analysis_2")
#
#     figure_settings_on_analysis_verticalLayout = QtWidgets.QVBoxLayout(verticalLayoutWidget_on_analysis_2)
#     figure_settings_on_analysis_verticalLayout.setContentsMargins(0, 0, 0, 0)
#     figure_settings_on_analysis_verticalLayout.setObjectName("figure_settings_on_analysis_verticalLayout")
#
#     layout.addWidget(analysis_settings_groupBox, 1, 2, 5, 1)
#
#
#
#     pymol_width_label_on_analysis = QtWidgets.QLabel(verticalLayoutWidget_on_analysis_2)
#     pymol_width_label_on_analysis.setText("Width: 1200")
#     sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
#     sizePolicy.setHorizontalStretch(0)
#     sizePolicy.setVerticalStretch(0)
#     sizePolicy.setHeightForWidth(pymol_width_label_on_analysis.sizePolicy().hasHeightForWidth())
#     pymol_width_label_on_analysis.setSizePolicy(sizePolicy)
#     pymol_width_label_on_analysis.setMinimumSize(QtCore.QSize(50, 22))
#     pymol_width_label_on_analysis.setMaximumSize(QtCore.QSize(16777215, 22))
#     pymol_width_label_on_analysis.setStyleSheet("QLabel {"
#                                                 "    background-color: rgb(27, 29, 35);"
#                                                 "    border-radius: 5px;"
#                                                 "    border: 2px solid rgb(27, 29, 35);"
#                                                 "    padding: 1px 1px 1px 1px;"
#                                                 "    border-bottom-color: rgb(157, 90, 198);"
#                                                 "}"
#                                                 ""
#                                                 "QLabel:hover{"
#                                                 "    border: 2px solid rgb(64, 71, 88);"
#                                                 "    selection-color: rgb(127, 5, 64);"
#                                                 "}")
#     pymol_width_label_on_analysis.setObjectName("pymol_width_label_on_analysis")
#     figure_settings_on_analysis_verticalLayout.addWidget(pymol_width_label_on_analysis)
#
#     width_horizontalSlider_on_analysis = QtWidgets.QSlider(verticalLayoutWidget_on_analysis_2)
#     sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
#     sizePolicy.setHorizontalStretch(0)
#     sizePolicy.setVerticalStretch(0)
#     sizePolicy.setHeightForWidth(width_horizontalSlider_on_analysis.sizePolicy().hasHeightForWidth())
#     width_horizontalSlider_on_analysis.setSizePolicy(sizePolicy)
#     width_horizontalSlider_on_analysis.setStyleSheet(
#         "QSlider::handle:horizontal {"
#         "    background-color:  rgb(255, 17, 100);"
#         "    border: 2px solid;"
#         "    width: 8px;"
#         "    margin: -15px 0px;"
#         "}"
#         ""
#         "QSlider:horizontal {"
#         "    min-height: 20px;"
#         "}"
#         ""
#         "QSlider::groove:horizontal {"
#         "    height: 1px;"
#         "    background-color: rgb(110, 105, 225);"
#         "    border: 1px solid;"
#         "    height: 5px;"
#         "    margin: 0px;"
#         "    border-radius: 5px;"
#         "}"
#         ""
#         "QSlider::handle:horizontal {"
#         "    width: 10px;"
#         "    margin-top: -10px;"
#         "    margin-bottom: -10px;"
#         "    border-radius: 5px;"
#         "    background-color: rgb(255, 17, 100);"
#         "    border: 2px solid;"
#         "}"
#         ""
#         "QSlider::handle:horizontal:hover {"
#         "    background-color: rgb(22, 200, 244);"
#         "}"
#     )
#     width_horizontalSlider_on_analysis.setMinimum(800)
#     width_horizontalSlider_on_analysis.setMaximum(1920)
#     width_horizontalSlider_on_analysis.setSingleStep(20)
#     width_horizontalSlider_on_analysis.setPageStep(20)
#     width_horizontalSlider_on_analysis.setProperty("value", 1200)
#     width_horizontalSlider_on_analysis.setOrientation(QtCore.Qt.Horizontal)
#     width_horizontalSlider_on_analysis.setInvertedAppearance(False)
#     width_horizontalSlider_on_analysis.setInvertedControls(False)
#     width_horizontalSlider_on_analysis.setTickPosition(QtWidgets.QSlider.NoTicks)
#     width_horizontalSlider_on_analysis.setTickInterval(20)
#     width_horizontalSlider_on_analysis.setObjectName("width_horizontalSlider_on_analysis")
#     figure_settings_on_analysis_verticalLayout.addWidget(width_horizontalSlider_on_analysis)
#
#     pymol_height_label_on_analysis = QtWidgets.QLabel(verticalLayoutWidget_on_analysis_2)
#     pymol_height_label_on_analysis.setText("Height: 1080")
#     sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
#     sizePolicy.setHorizontalStretch(0)
#     sizePolicy.setVerticalStretch(0)
#     sizePolicy.setHeightForWidth(pymol_height_label_on_analysis.sizePolicy().hasHeightForWidth())
#     pymol_height_label_on_analysis.setSizePolicy(sizePolicy)
#     pymol_height_label_on_analysis.setMinimumSize(QtCore.QSize(50, 22))
#     pymol_height_label_on_analysis.setMaximumSize(QtCore.QSize(16777215, 22))
#     pymol_height_label_on_analysis.setStyleSheet("QLabel {"
#                                                  "    background-color: rgb(27, 29, 35);"
#                                                  "    border-radius: 5px;"
#                                                  "    border: 2px solid rgb(27, 29, 35);"
#                                                  "    padding: 1px 1px 1px 1px;"
#                                                  "    border-bottom-color: rgb(157, 90, 198);"
#                                                  "}"
#                                                  ""
#                                                  "QLabel:hover{"
#                                                  "    border: 2px solid rgb(64, 71, 88);"
#                                                  "    selection-color: rgb(127, 5, 64);"
#                                                  "}")
#     pymol_height_label_on_analysis.setObjectName("pymol_height_label_on_analysis")
#     figure_settings_on_analysis_verticalLayout.addWidget(pymol_height_label_on_analysis)
#
#     height_horizontalSlider_on_analysis = QtWidgets.QSlider(verticalLayoutWidget_on_analysis_2)
#     sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
#     sizePolicy.setHorizontalStretch(0)
#     sizePolicy.setVerticalStretch(0)
#     sizePolicy.setHeightForWidth(height_horizontalSlider_on_analysis.sizePolicy().hasHeightForWidth())
#     height_horizontalSlider_on_analysis.setSizePolicy(sizePolicy)
#     height_horizontalSlider_on_analysis.setStyleSheet(
#         "QSlider::handle:horizontal {"
#         "    background-color:  rgb(255, 17, 100);"
#         "    border: 2px solid;"
#         "    width: 8px;"
#         "    margin: -15px 0px;"
#         "}"
#         ""
#         "QSlider:horizontal {"
#         "    min-height: 20px;"
#         "}"
#         ""
#         "QSlider::groove:horizontal {"
#         "    height: 1px;"
#         "    background-color: rgb(110, 105, 225);"
#         "    border: 1px solid;"
#         "    height: 5px;"
#         "    margin: 0px;"
#         "    border-radius: 5px;"
#         "}"
#         ""
#         "QSlider::handle:horizontal {"
#         "    width: 10px;"
#         "    margin-top: -10px;"
#         "    margin-bottom: -10px;"
#         "    border-radius: 5px;"
#         "    background-color: rgb(255, 17, 100);"
#         "    border: 2px solid;"
#         "}"
#         ""
#         "QSlider::handle:horizontal:hover {"
#         "    background-color: rgb(22, 200, 244);"
#         "}"
#     )
#     height_horizontalSlider_on_analysis.setMinimum(600)
#     height_horizontalSlider_on_analysis.setMaximum(1440)
#     height_horizontalSlider_on_analysis.setSingleStep(20)
#     height_horizontalSlider_on_analysis.setPageStep(20)
#     height_horizontalSlider_on_analysis.setProperty("value", 1080)
#     height_horizontalSlider_on_analysis.setOrientation(QtCore.Qt.Horizontal)
#     height_horizontalSlider_on_analysis.setInvertedAppearance(False)
#     height_horizontalSlider_on_analysis.setInvertedControls(False)
#     height_horizontalSlider_on_analysis.setTickPosition(QtWidgets.QSlider.NoTicks)
#     height_horizontalSlider_on_analysis.setTickInterval(20)
#     height_horizontalSlider_on_analysis.setObjectName("height_horizontalSlider_on_analysis")
#     figure_settings_on_analysis_verticalLayout.addWidget(height_horizontalSlider_on_analysis)
#
#     layout.addWidget(analysis_settings_groupBox, 1, 2, 5, 1)

