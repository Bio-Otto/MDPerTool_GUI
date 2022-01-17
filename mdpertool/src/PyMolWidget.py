# from pylab import *
# from PySide2.QtOpenGL import *
import os.path
from PySide2 import QtCore, QtWidgets
from PySide2.QtCore import Qt
from OpenGL.GL import *
import pymol2
from PySide2.QtOpenGL import QGLFormat, QGLWidget
from pymol.cgo import *
from pymol.vfont import plain
from pymol.cgo import *
from pymol.vfont import plain
from pymol import cmd, cgo, CmdException, selector, stored
from chempy import cpv

os.environ['QT_API'] = 'pyside2'
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
        # self._pymol.cmd.fragment("his")

        self.mol_name = None

    def initializeGL(self):
        """
        Reimplemented from QGLWidget
        Instance PyMOL _only_ when we're sure there's an OGL context up and running
        (i.e. in this method :-)
        """
        print(path)
        if not self._enableUi:
            self._pymol.cmd.set("internal_gui", 0)
            self._pymol.cmd.set("internal_feedback", 0)
            self._pymol.cmd.button("double_left", "None", "None")
            self._pymol.cmd.button("single_right", "None", "None")

        self._pymol.reshape(self.width(), self.height())
        self.resizeGL(self.width(), self.height())
        self._pymolProcess()

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

    def set_ss_figure(self):
        self._pymol.cmd.set('ray_opaque_background', 0)
        self._pymol.cmd.color('red', 'ss h')
        self._pymol.cmd.color('yellow', 'ss s')
        self._pymol.cmd.set('ray_trace_mode', 3)

    def get_png_figure(self, figure_name, width=1200, height=1200, dpi=150, ray=1):
        self._pymol.cmd.png(figure_name, width=width, height=height, dpi=dpi, ray=ray)

    def loadMolFile(self, mol_file):
        try:
            self._pymol.cmd.load(str(mol_file))
            self.mol_name = os.path.basename(mol_file).split('.')[0]
        except Exception as loadingError:
            print("An error occurred while loading the structure file :(\n", loadingError)

    def initial_pymol_visual(self):
        cgo = []
        axes = [[2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 2.0]]
        pos = [1.0, 1.0, -10.0]
        cyl_text(cgo, plain, pos, 'MDPerTool', 0.20, [1.0, 0.1, 0.1], axes=axes)
        pos = [0.0, -2.0, -15.0]
        wire_text(cgo, plain, pos, 'Version v0.1', axes)
        self._pymol.cmd.bg_color('0x2C313C')
        self._pymol.cmd.set("cgo_line_radius", 0.05)
        self._pymol.cmd.load_cgo(cgo, 'txt')
        self._pymol.cmd.zoom()

        try:
            self.loadMolFile(demo_pdb_path)
            self.mol_name = os.path.basename(demo_pdb_path).split('.')[0]
        except:
            self._pymol.cmd.fragment("phe")

        self._pymol.cmd.center()
        self._pymol.cmd.zoom()

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

        except Exception as expression:
            print(expression)

    def resi_label_add(self, selection):
        try:
            self._pymol.cmd.select(selection)
            # self._pymol.cmd.do('color red, ' + selection)

            label_selection = '%s and name ca' % selection
            self._pymol.cmd.set('label_color', 'green', label_selection)
            self._pymol.cmd.label(label_selection, '"%s-%s" % (resn, resi)')
            self._pymol.cmd.set('label_position', '(1, 1, 6)')
            self._pymol.cmd.set('label_size', '14')

        except Exception as expression:
            print(expression)

    def clear_all_labels(self):
        self._pymol.cmd.label('all', '')
        self._pymol.cmd.do('color green')

    def activate_navigation_tool(self):
        self._pymol.cmd.set("internal_gui", 1)
        self._pymol.cmd.set("internal_feedback", 1)
        self._pymol.cmd.button("double_left", "None", "None")
        self._pymol.cmd.button("single_right", "None", "None")
        self._pymol.reshape(self.width(), self.height())
        self.resizeGL(self.width(), self.height())
        self._pymolProcess()

    def deactivate_navigation_tool(self):
        self._pymol.cmd.set("internal_gui", 0)
        self._pymol.cmd.set("internal_feedback", 0)
        self._pymol.cmd.button("double_left", "None", "None")
        self._pymol.cmd.button("single_right", "None", "None")
        self._pymol.reshape(self.width(), self.height())
        self.resizeGL(self.width(), self.height())
        self._pymolProcess()

    # def show_energy_dissipation(self, response_time_file_path, spectrum='white_blue', mol=None,
    #                             display_color_ramp=True):
    #
    #     """
    #     :param self:
    #     :param response_time_file_path: Response Time csv file that include residue responses
    #     :param mol: Structure file path that will be display on
    #     :param display_color_ramp: Color bar will appear on screen
    #     :param spectrum: Transition colors in the structure
    #     :return:
    #     """
    #
    #     # load the protein
    #     print(self.mol_name)
    #     if mol is None:
    #         mol = self.mol_name
    #
    #     else:
    #         self.mol_name = mol
    #     # open the file of new values (just 1 column of numbers, one for each alpha carbon)
    #     inFile = open(response_time_file_path, 'r')
    #
    #     # create the global, stored array
    #     stored = []
    #
    #     # read the new B factors from file
    #     for line in inFile.readlines():
    #         stored.append(float(line))
    #
    #     # close the input file
    #     inFile.close()
    #
    #     # clear out the old B Factors
    #     self._pymol.cmd.alter('all', 'b=0.0')
    #
    #     # ---> obj = self._pymol.cmd.get_object_list(mol)[0]
    #
    #     counter = 0
    #     for line in stored:
    #         bfact = float(line)
    #
    #         self._pymol.cmd.alter("%s and resi %s and n. CA" % (str(self.mol_name), counter), "b=%s" % bfact)
    #         counter = counter + 1
    #
    #     self._pymol.cmd.spectrum('b', spectrum, minimum=0, maximum=len(stored))
    #     self._pymol.cmd.recolor()
    #
    #     if display_color_ramp:
    #         self._pymol.cmd.ramp_new("count", self.mol_name, [0, len(stored)], color=[[1.0, 1.0, 1.0], 'blue'])

    # def create_directed_arrows(self, atom1='pk1', atom2='pk2', radius=0.5, gap=0.0, hlength=-1, hradius=-1,
    #                            color='red', name=''):
    #
    #     """
    #         cgo_arrow [ atom1 [, atom2 [, radius [, gap [, hlength [, hradius [, color [, name ]]]]]]]]
    #         cgo_arrow A/PNK`301/C13, A/PNK`301/C7, gap=1.0, radius=0.15, hradius=0.5, hlength=0.5
    #         cgo_arrow A/TYR`7/C, A/ASN`2/N, color=red red, hradius=1.0, hlength=1.0, gap=1.0, radius=0.35
    #         atom1 = string: single atom selection or list of 3 floats {default: pk1}
    #         atom2 = string: single atom selection or list of 3 floats {default: pk2}
    #         radius = float: arrow radius {default: 0.5}
    #         gap = float: gap between arrow tips and the two atoms {default: 0.0}
    #         hlength = float: length of head
    #         hradius = float: radius of head
    #         color = string: one or two color names {default: blue red}
    #         name = string: name of CGO object
    #     """
    #
    #     radius, gap = float(radius), float(gap)
    #     hlength, hradius = float(hlength), float(hradius)
    #
    #     try:
    #         color1, color2 = color.split()
    #     except:
    #         color1 = color2 = color
    #     color1 = list(self._pymol.cmd.get_color_tuple(color1))
    #     color2 = list(self._pymol.cmd.get_color_tuple(color2))
    #
    #     def get_coord(v):
    #         if not isinstance(v, str):
    #             return v
    #         if v.startswith('['):
    #             return cmd.safe_list_eval(v)
    #         return cmd.get_atom_coords(v)
    #
    #     xyz1 = get_coord(atom1)
    #     xyz2 = get_coord(atom2)
    #     normal = cpv.normalize(cpv.sub(xyz1, xyz2))
    #
    #     if hlength < 0:
    #         hlength = radius * 3.0
    #     if hradius < 0:
    #         hradius = hlength * 0.6
    #
    #     if gap:
    #         diff = cpv.scale(normal, gap)
    #         xyz1 = cpv.sub(xyz1, diff)
    #         xyz2 = cpv.add(xyz2, diff)
    #
    #     xyz3 = cpv.add(cpv.scale(normal, hlength), xyz2)
    #
    #     obj = [cgo.CYLINDER] + xyz1 + xyz3 + [radius] + color1 + color2 + \
    #           [cgo.CONE] + xyz3 + xyz2 + [hradius, 0.0] + color2 + color2 + \
    #           [1.0, 0.0]
    #
    #     if not name:
    #         name = self._pymol.cmd.get_unused_name('arrow')
    #
    #     self._pymol.cmd.load_cgo(obj, name)
    #     self._pymol.cmd.group("Arrows", name)

    def MinMaxScaler(self, X, min, max, top):
        return (X - min) / (max - min) * top

    def show_energy_dissipation(self, response_time_file_path, spectrum='red_white', mol=None,
                                display_color_ramp=True):

        """
            :param self:
            :param response_time_file_path: Response Time csv file that include residue responses
            :param mol: Structure file path that will be display on
            :param display_color_ramp: Color bar will appear on screen
            :param spectrum: Transition colors in the structure
            :return:
        """
        self._pymol.cmd.show('surface')
        self._pymol.cmd.set('surface_quality', 1)
        self._pymol.cmd.set('surface_color', "white")
        self._pymol.cmd.set('transparency', 0.55)

        self._pymol.stored.residues = []
        self._pymol.cmd.iterate('name ca', 'stored.residues.append(resi)', _self=self._pymol.cmd)

        start_res_number = [int(x) for x in self._pymol.stored.residues][0]

        # load the protein
        if mol is None:
            mol = self.mol_name

        else:
            self.mol_name = mol
        # open the file of new values (just 1 column of numbers, one for each alpha carbon)
        inFile = open(response_time_file_path, 'r')

        # create the global, stored array
        stored = []
        normalize_stored = []
        # read the new B factors from file
        for line in inFile.readlines():
            stored.append(float(line))
        # close the input file
        inFile.close()

        MIN = min(stored)
        MAX = max(stored)
        # Max Min Scaler function

        for i in stored:
            norm_val = round(self.MinMaxScaler(i, MIN, MAX, top=len(stored)), 2)
            normalize_stored.append(norm_val)

        stored = normalize_stored

        # clear out the old B Factors
        self._pymol.cmd.alter('all', 'b=0.0')

        # ---> obj = self._pymol.cmd.get_object_list(mol)[0]
        counter = start_res_number
        for line in stored:
            bfact = float(line)
            self._pymol.cmd.alter("%s and resi %s and n. CA" % (str(self.mol_name), counter), "b=%s" % bfact)
            counter = counter + 1

        self._pymol.cmd.spectrum('b', spectrum, minimum=0, maximum=max(stored))
        self._pymol.cmd.recolor()

        if display_color_ramp:
            self._pymol.cmd.ramp_new("count", self.mol_name, [0, len(stored)], color=['red', 'white'])

        # self._pymol.cmd.save("session.pse")
        # self._pymol.cmd.save('aa.pdb')

    def create_directed_arrows(self, atom1='pk1', atom2='pk2', radius=0.5, gap=0.0, hlength=-1, hradius=-1,
                               color='red', name='', shortest_path=False):

        """
            cgo_arrow [ atom1 [, atom2 [, radius [, gap [, hlength [, hradius [, color [, name ]]]]]]]]
            cgo_arrow A/PNK`301/C13, A/PNK`301/C7, gap=1.0, radius=0.15, hradius=0.5, hlength=0.5
            cgo_arrow A/TYR`7/C, A/ASN`2/N, color=red red, hradius=1.0, hlength=1.0, gap=1.0, radius=0.35
            atom1 = string: single atom selection or list of 3 floats {default: pk1}
            atom2 = string: single atom selection or list of 3 floats {default: pk2}
            radius = float: arrow radius {default: 0.5}
            gap = float: gap between arrow tips and the two atoms {default: 0.0}
            hlength = float: length of head
            hradius = float: radius of head
            color = string: one or two color names {default: blue red}
            name = string: name of CGO object
        """

        try:
            aa = self._pymol.cmd.get_names("nongroup_objects", enabled_only=1)
            clean_ls = [i for i in aa if i.startswith('arrow')]
            for i in clean_ls:
                cc = self._pymol.cmd.get_model('(%s)' % i, 1).get_coord_list()
                print(cc)

        except:
            pass

        if shortest_path:
            color = 'blue'

        else:
            pass

        radius, gap = float(radius), float(gap)
        hlength, hradius = float(hlength), float(hradius)

        try:
            color1, color2 = color.split()
        except:
            color1 = color2 = color
        color1 = list(self._pymol.cmd.get_color_tuple(color1))
        color2 = list(self._pymol.cmd.get_color_tuple(color2))

        def get_coord(v):
            if not isinstance(v, str):
                return v
            if v.startswith('['):
                return cmd.safe_list_eval(v)
            return cmd.get_atom_coords(v)

        xyz1 = get_coord(atom1)
        xyz2 = get_coord(atom2)
        normal = cpv.normalize(cpv.sub(xyz1, xyz2))

        if hlength < 0:
            hlength = radius * 3.0
        if hradius < 0:
            hradius = hlength * 0.6

        if gap:
            diff = cpv.scale(normal, gap)
            xyz1 = cpv.sub(xyz1, diff)
            xyz2 = cpv.add(xyz2, diff)

        xyz3 = cpv.add(cpv.scale(normal, hlength), xyz2)

        obj = [cgo.CYLINDER] + xyz1 + xyz3 + [radius] + color1 + color2 + \
              [cgo.CONE] + xyz3 + xyz2 + [hradius, 0.0] + color2 + color2 + \
              [1.0, 0.0]

        if not shortest_path:

            if name == '':
                new_name = self._pymol.cmd.get_unused_name('arrow')
                self._pymol.cmd.load_cgo(obj, new_name)
                self._pymol.cmd.group('Arrows', name)
                # self._pymol.cmd.load_cgo(obj, name)
                # self._pymol.cmd.group("Arrows", name)

            if name != '':
                new_name = self._pymol.cmd.get_unused_name(str(name))
                self._pymol.cmd.load_cgo(obj, new_name)
                self._pymol.cmd.group(name + 's', new_name)

        if shortest_path:
            hradius = -0.5
            radius = 0.10
            xyz1 = get_coord(atom1)
            xyz2 = get_coord(atom2)
            normal = cpv.normalize(cpv.sub(xyz1, xyz2))

            if gap:
                diff = cpv.scale(normal, gap)
                xyz1 = cpv.sub(xyz1, diff)
                xyz2 = cpv.add(xyz2, diff)

            xyz3 = cpv.add(cpv.scale(normal, hlength), xyz2)

            obj = [cgo.CYLINDER] + xyz1 + xyz3 + [radius] + color1 + color2 + \
                  [cgo.CONE] + xyz3 + xyz2 + [hradius, 0.0] + color2 + color2 + \
                  [1.0, 0.0]

            name = self._pymol.cmd.get_unused_name('path')
            self._pymol.cmd.load_cgo(obj, name)
            self._pymol.cmd.group("Sh_Paths", name)

        self.update()
        # cmd.extend('cgo_arrow', cgo_arrow)

    def create_interacting_Residues(self, atom1='pk1', atom2='pk2', radius=0.5, gap=0.0, hlength=-1, hradius=-1,
                                    color='red', name=''):

        try:
            aa = self._pymol.cmd.get_names("nongroup_objects", enabled_only=1)
            clean_ls = [i for i in aa if i.startswith('interact')]
            for i in clean_ls:
                self._pymol.cmd.get_model('(%s)' % i, 1).get_coord_list()

        except:
            pass

        radius, gap = float(radius), float(gap)
        hlength = float(hlength)

        try:
            color1, color2 = color.split()
        except:
            color1 = color2 = color
        color1 = list(self._pymol.cmd.get_color_tuple(color1))
        color2 = list(self._pymol.cmd.get_color_tuple(color2))

        def get_coord(v):
            if not isinstance(v, str):
                return v
            if v.startswith('['):
                return cmd.safe_list_eval(v)
            return cmd.get_atom_coords(v)

        xyz1 = get_coord(atom1)
        xyz2 = get_coord(atom2)
        normal = cpv.normalize(cpv.sub(xyz1, xyz2))

        if hlength < 0:
            hlength = radius * 3.0

        if gap:
            diff = cpv.scale(normal, gap)
            xyz1 = cpv.sub(xyz1, diff)
            xyz2 = cpv.add(xyz2, diff)

        xyz3 = cpv.add(cpv.scale(normal, hlength), xyz2)
        obj = [cgo.CYLINDER] + xyz1 + xyz3 + [radius] + color1 + color2

        if not name:
            name = self._pymol.cmd.get_unused_name('interact')

        self._pymol.cmd.load_cgo(obj, name)
        self._pymol.cmd.group("Interacts", name)

        self.update()

    def show_ligand_polar_interactions(self):
        self._pymol.preset.ligands(selection='all', _self=self._pymol.cmd)
        self._pymol.cmd.show("labels", 'chin_ins_pose_1_pol_conts')

        # objs = self._pymol.cmd.get_object_list('all')
        # objs2 = self._pymol.cmd.get_names('objects', 0, selection='(all)')
        # print(objs2)

        # cmd.extend('cgo_arrow', cgo_arrow)
# if __name__ == "__main__":
#     app = QtWidgets.QApplication()
#     window = PymolQtWidget()
#     window.show()
#     app.exec_()
#     window = PymolQtWidget()
#     window.show()
