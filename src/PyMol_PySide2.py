from PySide2 import QtCore, QtWidgets
from PySide2.QtCore import Qt
from PySide2.QtOpenGL import QGLWidget, QGLFormat
import os
import pymol2

os.environ['QT_API'] = 'pyside2'

buttonMap = {
    Qt.LeftButton: 0,
    Qt.MidButton: 1,
    Qt.RightButton: 2,
}


class PymolQtWidget(QGLWidget):

    def __init__(self, parent=None):
        self._enableUi = False
        f = QGLFormat()
        f.setStencil(True)
        f.setRgba(True)
        f.setDepth(True)
        f.setDoubleBuffer(True)

        self._pymol = pymol2.PyMOL()
        self._pymol.start()
        super(PymolQtWidget, self).__init__(f, parent)

        self._timer = QtCore.QTimer()
        self._timer.setSingleShot(True)
        self._timer.timeout.connect(self._pymolProcess)

        # Load a demo molecule
        self._pymol.cmd.fragment("his")

    def initializeGL(self):
        """
        Reimplemented from QGLWidget
        Instance PyMOL _only_ when we're sure there's an OGL context up and running
        (i.e. in this method :-)
        """
        if not self._enableUi:
            self._pymol.cmd.set("internal_gui", 0)
            self._pymol.cmd.set("internal_feedback", 0)
            self._pymol.cmd.button("double_left", "None", "None")
            self._pymol.cmd.button("single_right", "None", "None")

        self._pymol.reshape(self.width(), self.height())
        self.resizeGL(self.width(), self.height())
        self._pymolProcess()

    def paintGL(self):
        self._doIdle()
        self._pymol.draw()

    def reinitialize(self):
        self._pymol.cmd.reinitialize()

    def resizeGL(self, w, h):
        self._pymol.reshape(w, h, True)
        self._pymolProcess()

    def _pymolProcess(self):
        self._doIdle()
        self.update()

    def mouseMoveEvent(self, ev):
        self._pymol.drag(ev.x(), self.height() - ev.y(), 0)
        self._pymolProcess()

    def mousePressEvent(self, ev):
        if not self._enableUi:
            self._pymol.cmd.button("double_left", "None", "None")
            self._pymol.cmd.button("single_right", "None", "None")
        self._pymol.button(buttonMap[ev.button()], 0, ev.x(), self.height() - ev.y(), 0)
        self._pymolProcess()

    def mouseReleaseEvent(self, ev):
        self._pymol.button(buttonMap[ev.button()], 1, ev.x(), self.height() - ev.y(), 0)
        self._pymolProcess()

    def wheelEvent(self, ev):
        self.button = 3
        self._pymol.button(self.button, 0, ev.x(), ev.y(), 0)
        self._pymolProcess()
        self._timer.start(0)

    def _doIdle(self):
        if self._pymol.idle():
            self._timer.start(0)


# if __name__ == "__main__":
#     app = QtWidgets.QApplication()
#     window = PymolQtWidget()
#     window.show()
#     app.exec_()
