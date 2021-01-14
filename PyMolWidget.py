from pylab import *
from PyQt5.QtOpenGL import *
from PyQt5 import QtCore, QtWidgets
from PyQt5.Qt import Qt
from OpenGL.GL import *
import pymol2
from PyQt5.QtOpenGL import QGLFormat
from pymol.cgo import *
from pymol.vfont import plain

buttonMap = {
    Qt.LeftButton: 0,
    Qt.MidButton: 1,
    Qt.RightButton: 2,
}
path = os.getcwd()
print(path)
demo_pdb_path = path + '/Download/1aki.pdb'


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

    def initializeGL(self):
        """
        Reimplemented from QGLWidget

        Instance PyMOL _only_ when we're sure there's an OGL context up and running
        (i.e. in this method :-)
        """

        # self._pymol.start()

        if not self._enableUi:
            self._pymol.cmd.set("internal_gui", 0)
            self._pymol.cmd.set("internal_feedback", 0)
            self._pymol.cmd.button("double_left", "None", "None")
            self._pymol.cmd.button("single_right", "None", "None")

        self._pymol.reshape(self.width(), self.height())
        self.resizeGL(self.width(), self.height())
        self._pymolProcess()
        # self._pymol.start()

    def paintGL(self):
        glViewport(0, 0, self.width(), self.height())
        self._doIdle()
        # self._pymol.idle()
        self._pymol.draw()

    def reinitialize(self):
        self._pymol.cmd.reinitialize()

    def resizeGL(self, w, h):
        self._pymol.reshape(w, h, True)
        self._pymolProcess()

    def loadMolFile(self, mol_file):
        print("geldi")

        self._pymol.cmd.load(str(mol_file))
        # self._pymol.cmd.bg_color('grey')

    def initial_pymol_visual(self):
        cgo = []
        axes = [[2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 2.0]]
        pos = [0.0, 0.0, 0.0]
        # wire_text(cgo, plain, pos, 'MDPerTool', axes)
        cyl_text(cgo, plain, pos, 'MDPerTool', 0.20, [1.0, 0.1, 0.1], axes=axes)

        pos = [0.0, -3.0, 0.0]
        wire_text(cgo, plain, pos, 'Version v0.1', axes)
        self._pymol.cmd.bg_color('0x2C313C')
        self._pymol.cmd.set("cgo_line_radius", 0.05)
        self._pymol.cmd.load_cgo(cgo, 'txt')
        self._pymol.cmd.zoom("all", -10.0)
        self.loadMolFile(demo_pdb_path)
        self._pymol.cmd.center("all", 0)
        # self._pymol.cmd.cgo.COLOR([0.2,0.3,0.4])
        # self._pymol.cmd.cgo.color(0.2, 0.3, 0.4)
        # self._pymol.cmd.full_screen('on')

    def selection_color(self, selection):

            # self._pymol.cmd.color('red', 'resi %s' %s[3])
        print(selection[3:-1])
        self.resicolor('resi %s' % selection[3:-1])


    def resicolor(self, selection):

        try:
            # self._pymol.cmd.full_screen('on')

            self._pymol.cmd.do('color green')
            self._pymol.cmd.select(selection)
            self._pymol.cmd.do('color red, ' + selection)
            self._pymol.cmd.hide('all')
            self._pymol.cmd.show('cartoon')
            label_selection = '%s and name ca' % selection
            self._pymol.cmd.set('label_color', 'blue', label_selection)
            self._pymol.cmd.label(label_selection, '"%s-%s" % (resn, resi)')
            self._pymol.cmd.set('label_position', '(3, 2, 1)')
            self._pymol.cmd.set('label_size', '20')
        except:
            print("problemm")

    def clear_all_labels(self):
        self._pymol.cmd.label('all', '')

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
        # self.button = 3 if ev.delta() > 0 else 4
        # self._pymol.button(button, 0, ev.x(), ev.y(), 0)
        # self._pymolProcess()
        self.button = 3
        self._pymol.button(self.button, 0, ev.x(), ev.y(), 0)
        self._pymolProcess()
        self._timer.start(0)

    def _doIdle(self):
        if self._pymol.idle():
            self._timer.start(0)

    def set_ss_figure(self):
        self._pymol.cmd.set('ray_opaque_background', 0)
        self._pymol.cmd.color('red', 'ss h')
        self._pymol.cmd.color('yellow', 'ss s')
        self._pymol.cmd.set('ray_trace_mode', 3)

    def get_png_figure(self):
        # self._pymol.cmd.full_screen('on')
        self._pymol.cmd.png("figure.png", width=1200, height=1200, dpi=300, ray=1)
