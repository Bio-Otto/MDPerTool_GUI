from PySide2 import QtCore, QtWidgets, QtWebEngineWidgets
import os
#from pyvis import network as net
from io import StringIO


class VisJS_QtWidget(QtWidgets.QWidget):
    def __init__(self, parent=None):
        super(VisJS_QtWidget, self).__init__(parent)

        self.html_file = None
        self.network = None
        self.m_output = QtWebEngineWidgets.QWebEngineView()

        # layout = QtWidgets.QVBoxLayout(self)
        # layout.addWidget(self.m_output)
        # self.resize(640, 480)

    def load_network_component(self, network, html_file='2d_network.html'):
        self.network = network
        self.html_file = html_file

    def __call__(self):
        g = net.Network(height='100%', width='70%', directed=True,
                        notebook=True, bgcolor='#2c313c', font_color='white')

        g.from_nx(self.network)
        g.show_buttons(filter_=['physics'])

        g.show(self.html_file)
        self.m_output.load(QtCore.QUrl().fromLocalFile(os.path.abspath(self.html_file)))
        # self.m_output.show()

