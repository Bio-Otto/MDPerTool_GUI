from PyQt5 import QtWidgets, uic, QtCore, QtGui
from pyqtgraph import PlotWidget, plot
import pyqtgraph as pg
import sys  # We need sys so that we can pass argv to QApplication
import os
import numpy as np
from PyQt5.QtCore import QTimer, QDateTime


class MainWindow(QtWidgets.QMainWindow):
    global curve, data, p6

    def __init__(self, *args, **kwargs):
        super(MainWindow, self).__init__(*args, **kwargs)

        # Load the UI Page
        uic.loadUi('MAIN_GUI.ui', self)

        #     self.graphWidget.setTitle("Temperature / Time", color="b", size="20pt") # Graph Title
        #
        #     self.pen = pg.mkPen(color=(0, 255, 0))  # plot line color
        #
        #
        #
        #     ## AXIS LABELS
        #     self.graphWidget.setLabel('left', "<span style=\"color:red;font-size:20px\">Temperature (°C)</span>")
        #     self.graphWidget.setLabel('bottom', "<span style=\"color:red;font-size:20px\">Hour (H)</span>")
        #
        #     self.plot([1, 2, 3, 4, 5, 6, 7, 8, 9, 10], [30, 32, 34, 32, 33, 31, 29, 32, 35, 45], self.pen,
        #               symbol='+')  # symbol='+'
        #     self.graphWidget.setBackground((100, 50, 255, 25))  # RGBA (A = alpha opacity) (100, 50, 255, 25)
        #
        # def plot(self, hour, temperature, pen, symbol=None, symbolBrush='b', symbol_size=20):
        #     self.graphWidget.plot(hour, temperature, pen=pen, symbol=symbol, symbolBrush=symbolBrush, symbolSize=symbol_size)
        #     # self.graphWidget.setXRange(290,320)
        #     # self.graphWidget.setYRange(290,320)
        # -*- coding: utf-8 -*-
        """
        This example demonstrates many of the 2D plotting capabilities
        in pyqtgraph. All of the plots may be panned/scaled by dragging with 
        the left/right mouse buttons. Right click on any plot to show a context menu.
        """

        win = pg.GraphicsLayoutWidget(show=True, title="Basic plotting examples")
        v_box = self.verticalLayout_16.addWidget(win)
        # win.resize(1000, 600)
        win.setWindowTitle('pyqtgraph example: Plotting')

        # Enable antialiasing for prettier plots
        pg.setConfigOptions(antialias=True)

        p1 = win.addPlot(title="Basic array plotting", y=np.random.normal(size=100))

        p2 = win.addPlot(title="Multiple curves")
        p2.plot(np.random.normal(size=100), pen=(255, 0, 0), name="Red curve")
        p2.plot(np.random.normal(size=110) + 5, pen=(0, 255, 0), name="Green curve")
        p2.plot(np.random.normal(size=120) + 10, pen=(0, 0, 255), name="Blue curve")

        p3 = win.addPlot(title="Drawing with points")
        p3.plot(np.random.normal(size=100), pen=(200, 200, 200), symbolBrush=(255, 0, 0), symbolPen='w')

        win.nextRow()

        p4 = win.addPlot(title="Parametric, grid enabled")
        x = np.cos(np.linspace(0, 2 * np.pi, 1000))
        y = np.sin(np.linspace(0, 4 * np.pi, 1000))
        p4.plot(x, y)
        p4.showGrid(x=True, y=True)

        p5 = win.addPlot(title="Scatter plot, axis labels, log scale")
        x = np.random.normal(size=1000) * 1e-5
        y = x * 1000 + 0.005 * np.random.normal(size=1000)
        y -= y.min() - 1.0
        mask = x > 1e-15
        x = x[mask]
        y = y[mask]
        p5.plot(x, y, pen=None, symbol='t', symbolPen=None, symbolSize=10, symbolBrush=(100, 100, 255, 50))
        p5.setLabel('left', "Y Axis", units='A')
        p5.setLabel('bottom', "Y Axis", units='s')
        p5.setLogMode(x=True, y=False)

        p6 = win.addPlot(title="Updating plot")
        curve = p6.plot(pen='y')
        data = np.random.normal(size=(300, 320))
        ptr = 0

        def update():
            nonlocal ptr
            print("geldi")

            curve.setData(data[ptr % 10])
            if ptr == 0:
                p6.enableAutoRange('xy', False)  ## stop auto-scaling after the first data set is plotted
            ptr += 1

        # timer = QtCore.QTimer()
        # timer.timeout.connect(update)
        # timer.start(50)

        self.timer = QTimer(self)
        self.timer.timeout.connect(update)
        self.timer.start(50)

        win.nextRow()

        p7 = win.addPlot(title="Filled plot, axis disabled")
        y = np.sin(np.linspace(0, 10, 1000)) + np.random.normal(size=1000, scale=0.1)
        p7.plot(y, fillLevel=-0.3, brush=(50, 50, 200, 100))
        p7.showAxis('bottom', False)

        x2 = np.linspace(-100, 100, 1000)
        data2 = np.sin(x2) / x2
        p8 = win.addPlot(title="Region Selection")
        p8.plot(data2, pen=(255, 255, 255, 200))
        lr = pg.LinearRegionItem([400, 700])
        lr.setZValue(-10)
        p8.addItem(lr)

        p9 = win.addPlot(title="Zoom on selected region")
        p9.plot(data2)

        def updatePlot():
            p9.setXRange(*lr.getRegion(), padding=0)

        def updateRegion():
            lr.setRegion(p9.getViewBox().viewRange()[0])

        lr.sigRegionChanged.connect(updatePlot)
        p9.sigXRangeChanged.connect(updateRegion)
        updatePlot()


# def main():
#     app = QtWidgets.QApplication(sys.argv)
#     main = MainWindow()
#     main.show()
#     sys.exit(app.exec_())
#
#
# if __name__ == '__main__':
#     main()

from PyQt5.QtWidgets import QPushButton, QApplication, QVBoxLayout, QWidget, QGroupBox, QHBoxLayout, QListWidget, \
    QFileDialog, QSizePolicy, QSpacerItem


class MyApp(object):
    def __init__(self):
        super(MyApp, self).__init__()
        self.mainWidget = QWidget()
        self.mainLayout = QVBoxLayout()
        self.mainWidget.setLayout(self.mainLayout)

        self.hLayout = QHBoxLayout()
        self.mainLayout.insertLayout(0, self.hLayout)
        self.listA = QListWidget()
        for i in range(3):
            self.listA.addItem('Item ' + str(i))
        self.hLayout.addWidget(self.listA)

        self.buttonGroupbox = QGroupBox()
        self.buttonlayout = QVBoxLayout()
        self.buttonGroupbox.setLayout(self.buttonlayout)

        okButton = QPushButton('Remove Selected')
        okButton.clicked.connect(self.removeSel)
        self.buttonlayout.addWidget(okButton)

        self.mainLayout.addWidget(self.buttonGroupbox)
        self.mainWidget.show()
        sys.exit(app.exec_())

    def removeSel(self):
        listItems = self.listA.selectedItems()
        if not listItems: return
        for item in listItems:
            self.listA.takeItem(self.listA.row(item))

    def add_item(self, res_num):
        self.listA.addItem(str(res_num))


# if __name__ == '__main__':
#     app = QApplication(sys.argv)
#     MyApp()


# This is the example code for loading files and content inside the file to QtGui.QListWidget
# It is PyQt4, but you can try with PyQt5 with small changes.
# If your are not expecting this answer, sorry.


class Window(QWidget):
    def __init__(self, parent=None):

        super(Window, self).__init__(parent)

        self.verticalLayout = QVBoxLayout(self)
        self.verticalLayout.setObjectName('verticalLayout')

        self.horizontalLayout = QHBoxLayout()
        self.horizontalLayout.setObjectName('horizontalLayout')

        self.listWidget = QListWidget(self)
        self.listWidget.setObjectName('listView')
        self.listWidget.setAlternatingRowColors(True)
        self.horizontalLayout.addWidget(self.listWidget)

        self.verticalLayout1 = QVBoxLayout()
        self.verticalLayout1.setSpacing(10)
        self.verticalLayout1.setObjectName('verticalLayout')

        self.pushButton = QPushButton(self)
        self.pushButton.setObjectName('pushButton')
        self.pushButton.setText('Load File Content')

        self.pushButton_2 = QPushButton(self)
        self.pushButton_2.setObjectName('pushButton_2')
        self.pushButton_2.setText('Load File')

        self.verticalLayout1.addWidget(self.pushButton)
        self.verticalLayout1.addWidget(self.pushButton_2)
        spacerItem = QSpacerItem(20, 40, QSizePolicy.Minimum, QSizePolicy.Expanding)
        self.verticalLayout1.addItem(spacerItem)

        self.horizontalLayout.addLayout(self.verticalLayout1)
        self.verticalLayout.addLayout(self.horizontalLayout)

        self.pushButton.clicked.connect(self.loadFileContent)
        self.pushButton_2.clicked.connect(self.loadFiles)

    def loadFileContent(self):
        openFiles = QFileDialog.getOpenFileName(self, 'Open File', 'c:/', 'txt (*.txt)')
        if openFiles:
            data = open(str(openFiles), 'r')
            dataList = data.readlines()
            self.listWidget.clear()

            for eachLine in dataList:
                if len(eachLine.strip()) != 0:
                    self.listWidget.addItem(eachLine.strip())

    def loadFiles(self):
        getDirectory = QFileDialog.getExistingDirectory(self, 'Browse', 'C:/')

        if getDirectory:
            fileList = os.listdir(str(getDirectory))

            if fileList:
                self.listWidget.clear()
                for eachFile in fileList:
                    self.listWidget.addItem(eachFile)


''' if __name__ == '__main__':
    app = QApplication(sys.argv)
    w = Window()
    w.show()
    sys.exit(app.exec_())
 '''