import os.path
import sys
import threading
import time
from PySide2 import QtCore, QtWidgets
from PySide2.QtCore import Qt, QThread, Signal, QTimer, QCoreApplication
from OpenGL.GL import *
import pymol2
from PySide2.QtOpenGL import QGLFormat, QGLWidget
from pymol.cgo import *
from pymol.vfont import plain
from pymol import cmd, cgo, CmdException, selector, stored
from chempy import cpv
from src.pdb_intro_1aki import pdb_id_1aki
from lxml import etree
import xml.etree.ElementTree as ET
import numpy as np


class MonitorWorker(QThread):
    file_updated = Signal(str)

    def __init__(self, pert_file_path, response_threshold, last_mod_time):
        super().__init__()
        self.pert_file_path = pert_file_path
        self.response_threshold = response_threshold
        self.last_mod_time = last_mod_time
        self._stop_thread_flag = False

    def run(self):
        """Continuously monitors the file and emits a signal if it is updated."""
        while not self._stop_thread_flag:
            try:
                mod_time = os.path.getmtime(self.pert_file_path)
                if mod_time > self.last_mod_time:
                    self.last_mod_time = mod_time
                    self.file_updated.emit(self.pert_file_path)
                time.sleep(2)
            except Exception as e:
                print(f"[ERROR] Exception in MonitorWorker: {e}")
                break

    def stop(self):
        """Stops the thread and waits for it to finish."""
        self._stop_thread_flag = True
        self.wait()  # Wait for the thread to finish


    def stop(self):
        """Stops the thread and waits for it to finish."""
        self._stop_thread_flag = True
        self.wait()  # Wait for the thread to finish


class FileProcessorWorker(QThread):
    update_complete = Signal()

    def __init__(self, pymol_widget, pert_file_path, response_threshold):
        super().__init__()
        self.pymol_widget = pymol_widget
        self.pert_file_path = pert_file_path
        self.response_threshold = response_threshold

    def run(self):
        """Processes the file and emits a signal when done."""
        try:
            self.pymol_widget.process_file(self.pert_file_path, threshold=self.response_threshold)
            self.update_complete.emit()
        except Exception as e:
            print(f"[ERROR] Exception in FileProcessorWorker: {e}")


os.environ['QT_API'] = 'pyside2'


# Mapping PySide2 mouse buttons to PyMOL expected values
# Correctly map mouse buttons to PyMOL expected values
buttonMap = {
    Qt.LeftButton: 0,
    Qt.RightButton: 2,
    Qt.MidButton: 1,
}

path = os.getcwd()
demo_pdb_path = os.path.join(path, 'Download', '1aki.pdb')


class PymolQtWidget(QGLWidget):
    def __init__(self, parent=None):
        self._initialize_pymol_instance()
        super().__init__(self._get_gl_format(), parent)

        self._setup_timers()
        self._initialize_widget_properties()

    def _initialize_pymol_instance(self):
        """Initialize the PyMOL instance."""
        self._pymol = pymol2.PyMOL()
        self._pymol.start()

    def _get_gl_format(self):
        """Configure and return the OpenGL format."""
        f = QGLFormat()
        f.setStencil(True)
        f.setRgba(True)
        f.setDepth(True)
        f.setDoubleBuffer(True)
        return f

    def _setup_timers(self):
        """Setup timers for the widget."""
        self._timer = QTimer()
        self._timer.setSingleShot(True)
        self._timer.timeout.connect(self._pymolProcess)

    def _initialize_widget_properties(self):
        """Initialize widget properties."""
        self._enableUi = False
        self.cumulative_affected_atoms = set()
        self.ref_file_path = None
        self.total_atom_count = None
        self.last_mod_time = 0
        self.total_steps = 0
        self.colored_atoms = {}
        self.colors_set = {}
        self.mol_name = None
        self.monitor_thread = None
        self.processor_thread = None
        self.effected_atom_count = None

    def initializeGL(self):
        """Initialize OpenGL context and PyMOL settings."""
        if not self._enableUi:
            self._configure_pymol_ui()

        self._pymol.reshape(self.width(), self.height())
        self._pymolProcess()

    def _configure_pymol_ui(self):
        """Configure PyMOL UI settings."""
        self._pymol.cmd.set("internal_gui", 0)
        self._pymol.cmd.set("internal_feedback", 0)
        self._pymol.cmd.button("double_left", "None", "None")
        self._pymol.cmd.button("single_right", "None", "None")

    def paintGL(self):
        """Reimplement QGLWidget's paintGL to handle PyMOL drawing."""
        self._update_viewport()
        self._doIdle()
        self._pymol.draw()

    def _update_viewport(self):
        """Update the OpenGL viewport."""
        glViewport(0, 0, self.width(), self.height())

    def resizeGL(self, w, h):
        """Handle widget resizing."""
        self._pymol.reshape(w, h, True)
        self._pymolProcess()

    def _pymolProcess(self):
        """Process PyMOL updates."""
        self._doIdle()
        self.update()

    def _doIdle(self):
        """Keep PyMOL responsive by processing idle tasks."""
        if self._pymol.idle():
            self._timer.start(0)

    def reinitialize(self):
        """Reinitialize PyMOL."""
        self._pymol.cmd.reinitialize()

    def change_default_background(self):
        """Change the default background color in PyMOL."""
        self._pymol.cmd.bg_color('0x2C313C')

    def set_ss_figure(self):
        """Set up a secondary structure figure in PyMOL."""
        self._pymol.cmd.set('ray_opaque_background', 0)
        self._pymol.cmd.color('red', 'ss h')
        self._pymol.cmd.color('yellow', 'ss s')
        self._pymol.cmd.set('ray_trace_mode', 3)

    # Event Handlers
    def mouseMoveEvent(self, ev):
        """Handle mouse move events."""
        self._pymol.drag(ev.x(), self.height() - ev.y(), 0)
        self._pymolProcess()

    def mousePressEvent(self, ev):
        """Handle mouse press events."""
        if not self._enableUi:
            self._configure_pymol_ui()

        self._handle_mouse_event(ev, 0)

    def mouseReleaseEvent(self, ev):
        """Handle mouse release events."""
        self._handle_mouse_event(ev, 1)

    def _handle_mouse_event(self, ev, state):
        """Handle mouse events and update PyMOL accordingly."""
        button = ev.button()

        if button in buttonMap:
            pymol_button = buttonMap[button]
            self._pymol.button(pymol_button, state, ev.x(), self.height() - ev.y(), 0)
        else:
            return  # Ignore events with unmapped buttons

        self._pymolProcess()

    def wheelEvent(self, ev):
        """Handle mouse wheel events."""
        move_distance = 20 if ev.angleDelta().y() > 0 else -20  # Hareket mesafesi
        self._pymol.cmd.move('z', move_distance)  # Sahneyi Z ekseninde hareket ettir

        self._pymolProcess()  # PyMOL işlemini güncelleme

    def get_png_figure(self, figure_name, width=1200, height=1200, dpi=150, ray=1):
        self._pymol.cmd.png(figure_name, width=width, height=height, dpi=dpi, ray=ray)

    def loadMolFile(self, mol_file):
        try:
            self._pymol.cmd.load(str(mol_file))
            self.mol_name = os.path.splitext(os.path.basename(mol_file))[0]
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
        try:
            # self._pymol.cmd.color('red', 'resi %s' %s[3])
            self.resicolor('resi %s' % selection[3:-1])
        except Exception as ErrorExp:
            print("ERROR on PyMOL (selection_color): ", ErrorExp)

    def resicolor(self, selection):
        try:
            # self._pymol.cmd.full_screen('on')
            label_selection = '%s and name ca' % selection
            self._pymol.cmd.do('color green')
            self._pymol.cmd.select(selection)

            self._pymol.cmd.do('color red, ' + selection)
            self._pymol.cmd.hide('all')
            self._pymol.cmd.show('cartoon')

            self._pymol.cmd.set('label_color', 'blue', label_selection)
            # self._pymol.cmd.label("(%s)" % label_selection, '"%s-%s" % (resn, resi)')

            # self._pymol.cmd.set('label_position', '(3, 2, 1)')
            # self._pymol.cmd.set('label_size', '20')

        except Exception as expression:
            print("ERROR on PyMOL (resicolor): ", expression)

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
            print("ERROR on PyMOL (Labeling): ", expression)

    def clear_all_labels(self):
        self._pymol.cmd.label('all', '')
        self._pymol.cmd.do('color green')

    def activate_navigation_tool(self):
        self._pymol.cmd.set("internal_gui", 1)
        self._pymol.cmd.set("internal_feedback", 0)
        # self._pymol.cmd.button("double_left", "None", "None")
        # self._pymol.cmd.button("single_right", "None", "None")
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
        self._pymol.cmd.set('surface_color', "palecyan")
        self._pymol.cmd.set('transparency', 0.75)

        self._pymol.cmd.hide('surface')

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

        # if display_color_ramp:
        #     self._pymol.cmd.ramp_new("count", self.mol_name, [0, len(stored)], color=['red', 'white'])

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

        except:
            pass

        if shortest_path:
            color = 'magenta'

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
                self._pymol.cmd.disable(name + 's')  # This will disable (hide) the group

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

        radius, gap = float(radius), float(gap)
        hlength = float(hlength)

        try:
            color1, color2 = color.split()
        except ValueError:
            color1 = color2 = color
        color1 = list(self._pymol.cmd.get_color_tuple(color1))
        color2 = list(self._pymol.cmd.get_color_tuple(color2))

        def get_coord(v):
            try:
                if not isinstance(v, str):
                    return v
                if v.startswith('['):
                    return cmd.safe_list_eval(v)
                coords = cmd.get_atom_coords(v)
                if coords is None:
                    raise ValueError(f"Invalid selection or coordinates: {v}")
                return coords
            except Exception as e:
                raise ValueError(f"Error retrieving coordinates for {v}: {e}")

        try:
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
            self._pymol.cmd.disable("Interacts")  # This will disable (hide) the group
            self.update()
        except Exception as e:
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
            lineno = exc_tb.tb_lineno

    def show_ligand_polar_interactions(self):
        self._pymol.preset.ligands(selection='all', _self=self._pymol.cmd)
        self._pymol.cmd.show("labels", 'chin_ins_pose_1_pol_conts')

    def show_protein_in_cartoon_and_atoms(self, selection='all', color='blue'):
        """
        Display the protein in both cartoon and sticks representation in PyMOL with consistent colors for atoms and cartoon.

        Parameters:
        - selection: The selection string for the protein to be shown. Default is 'all' for all objects.
        - color: The color to be applied to both cartoon and sticks representations. Default is 'blue'.
        """
        try:
            # Show cartoon representation
            self._pymol.cmd.show('cartoon', selection)
            # Show sticks representation
            self._pymol.cmd.show('sticks', selection)

            # Set a consistent color for cartoon representation
            self._pymol.cmd.color(color, 'cartoon ' + selection)

            # Set the same color for sticks representation
            self._pymol.cmd.color(color, selection)

            # Rebuild the PyMOL view
            self._pymol.cmd.zoom(selection)
        except Exception as e:
            print("Error:", e)

    def parse_velocity_file(self, filename):
        tree = ET.parse(filename)
        root = tree.getroot()
        velocities = {}

        for velocity in root.findall('.//Velocity'):
            step = int(velocity.get('step', 0))
            atom_velocities = []
            for atom in velocity.findall('Atom'):
                x = float(atom.get('x', 0))
                y = float(atom.get('y', 0))
                z = float(atom.get('z', 0))
                atom_velocities.append([x, y, z])
            velocities[step] = np.array(atom_velocities)

        return velocities

    def find_affected_atoms(self, velocity_differences, threshold=0.05):
        """
        Identifies atoms whose velocity differences exceed a certain threshold.

        Parameters:
        - velocity_differences: A list of tuples where each tuple contains (step, atom_index, dx, dy, dz).
        - threshold: The magnitude threshold to consider an atom as affected.

        Returns:
        - affected_atoms: A dictionary where keys are step numbers and values are lists of affected atom indices.
        """
        affected_atoms = {}
        step_diffs = {}

        # Organize differences by step
        for step, atom_index, dx, dy, dz in velocity_differences:
            if step not in step_diffs:
                step_diffs[step] = []
            step_diffs[step].append((atom_index, dx, dy, dz))

        # Find affected atoms based on threshold
        for step, diffs in step_diffs.items():
            # Calculate magnitudes
            magnitudes = np.linalg.norm(np.array(diffs)[:, 1:], axis=1)  # Only take dx, dy, dz

            # Filter atoms based on threshold
            affected_indices = np.where(magnitudes > threshold)[0].tolist()
            affected_atoms[step] = [diffs[i][0] for i in affected_indices]

        return affected_atoms

    def get_residue_atom_indices(self):

        residue_atoms = {}
        # Modeli al
        model = self._pymol.cmd.get_model("all")

        for atom in model.atom:
            residue_name = f"{atom.resn}{atom.resi}"  # Residue ismi (resn: residue name, resi: residue index)
            if residue_name not in residue_atoms:
                residue_atoms[residue_name] = []
            residue_atoms[residue_name].append(atom.index)

        return residue_atoms

    def update_atom_colors(self, affected_atoms):

        self.total_steps += len(affected_atoms)

        for step, indices in affected_atoms.items():
            color_value = step / self.total_steps
            color = [1.0, 1.0 - color_value, 1.0 - color_value]
            color_name = f'step_{step}'
            self._pymol.cmd.set_color(color_name, color)

            self.cumulative_affected_atoms.update(indices)
            for atom_index in self.cumulative_affected_atoms:
                selection_name = f"index_{atom_index}"
                self._pymol.cmd.select(selection_name, f"index {atom_index}")
                self._pymol.cmd.color(color_name, selection_name)
                self._pymol.cmd.delete(selection_name)

    def calculate_velocity_differences(self, ref_velocities, pert_velocities):
        steps = ref_velocities.keys()
        velocity_differences = {}
        for step in steps:
            ref_data = ref_velocities.get(step, np.array([]))
            pert_data = pert_velocities.get(step, np.array([]))
            if ref_data.size > 0 and pert_data.size > 0:
                velocity_differences[step] = pert_data - ref_data
        return velocity_differences

    def get_last_steps(self, differences):
        """
        Get the total number of unique steps from the differences list.
        """
        unique_steps = set(step for step, _, _, _, _ in differences)
        # Benzersiz step sayısını geri döndürün
        return len(unique_steps)

    def monitor_live_file(self, pert_file_path, response_threshold=0.05):
        self.monitor_thread = MonitorWorker(pert_file_path, response_threshold, self.last_mod_time)
        self.monitor_thread.file_updated.connect(self.on_file_updated)
        self.show_protein_in_cartoon_and_atoms()
        self.monitor_thread.start()

    def on_file_updated(self, pert_file_path):

        if self.processor_thread is None or not self.processor_thread.isRunning():
            self.processor_thread = FileProcessorWorker(self, pert_file_path, self.monitor_thread.response_threshold)
            self.processor_thread.update_complete.connect(self._pymolProcess)
            self.processor_thread.start()

    def stop_monitoring(self):
        if self.monitor_thread is not None:
            self.monitor_thread.stop()

            self.cumulative_affected_atoms.clear()
            self.ref_file_path = None
            self.effected_atom_keeper = None
            self.last_mod_time = 0
            self.total_steps = 0
            self.colored_atoms.clear()
            self.colors_set.clear()
            self.total_atom_count = None
            self.effected_atom_count =None

    def process_file(self, pert_file_path, threshold=0.05):
        """
        Process the perturbed XML file and update PyMOL.
        """
        ref_velocities = self.parse_velocities(self.ref_file_path)
        pert_velocities = self.parse_velocities(pert_file_path)
        differences = self.compare_velocities(ref_velocities, pert_velocities)
        affected_atoms = self.find_affected_atoms(differences, threshold=threshold)
        self.update_colors(affected_atoms, max_step=len(ref_velocities))

        if self.total_atom_count is None:
            self.total_atom_count = len(ref_velocities)

    def parse_velocities(self, filename):
        """
        Parse velocities from an XML file and store them in a dictionary.

        Parameters:
        - filename (str): The path to the XML file containing velocity data.

        Returns:
        - dict: A dictionary where keys are step numbers and values are lists of tuples representing atom velocities (x, y, z).
        """

        velocities = {}

        try:
            # Iterate through the XML file and parse 'Velocity' elements
            for event, element in etree.iterparse(filename, events=('end',), tag='Velocity'):
                # Extract the step number from the element attributes
                step = int(element.get('step'))

                # Extract atom velocities from the 'Atom' sub-elements
                atoms = [
                    (float(atom.get('x')), float(atom.get('y')), float(atom.get('z')))
                    for atom in element.findall('Atom')
                ]

                # Store the velocities in the dictionary with step as the key
                velocities[step] = atoms

                # Clear the element to free memory
                element.clear()

        except etree.XMLSyntaxError:
            pass
            #print(f"XML Syntax Error: {e}")
            # Handle XML syntax errors if needed

        except Exception as e:
            print(f"Unexpected Error: {e}")

        # Initialize the total_atom_count if it is not already set
        if self.total_atom_count is None:
            # Use the number of atoms from the first step
            self.total_atom_count = len(atoms)

        return velocities

    def compare_velocities(self, ref_velocities, pert_velocities):
        """
        Compare the reference velocities with the perturbed velocities for each step.

        Parameters:
        - ref_velocities (dict): Dictionary where keys are step numbers and values are lists of reference atom velocities (x, y, z).
        - pert_velocities (dict): Dictionary where keys are step numbers and values are lists of perturbed atom velocities (x, y, z).

        Returns:
        - list: A list of tuples containing the step number, atom index, and differences in velocities (dx, dy, dz).
        """
        differences = []

        # Iterate through each step in the perturbed velocities
        for step, pert_atoms in pert_velocities.items():
            # Retrieve reference atoms for the current step; default to empty list if not present
            ref_atoms = ref_velocities.get(step, [])

            # Iterate through each perturbed atom
            for i, (px, py, pz) in enumerate(pert_atoms):
                # Get the reference atom's velocities; default to (0, 0, 0) if not present
                if i < len(ref_atoms):
                    rx, ry, rz = ref_atoms[i]
                else:
                    rx, ry, rz = 0, 0, 0

                # Calculate the differences in velocities
                dx = px - rx
                dy = py - ry
                dz = pz - rz

                # Record the step, atom index, and velocity differences
                differences.append((step, i, dx, dy, dz))

        return differences

    def get_color_for_step(self, step, max_step):
        """
        Generate a color for a given step based on its position relative to the maximum step.

        Parameters:
        - step (int): The current step number for which the color needs to be generated.
        - max_step (int): The maximum step number, which is used to normalize the step into a color gradient.

        Returns:
        - list: A list representing the RGB color values [red, green, blue].
                Returns None if the step is outside the valid range.
        """
        # Define the minimum step number
        min_step = 1

        # Check if the step is within the valid range
        if step < min_step or step > max_step:
            return None

        # Normalize the step to a range between 0 and 1
        normalized_step = (step - min_step) / (max_step - min_step)

        # Generate color based on the normalized step value
        # Starting with red, transitioning to white as the step number increases
        red = 1.0  # Red stays constant
        green = normalized_step  # Green increases with the step number
        blue = normalized_step  # Blue increases with the step number

        # Return the RGB color as a list
        return [red, green, blue]

    def update_colors(self, affected_atoms, max_step):
        """
        Update the colors of affected atoms in PyMOL based on their step number.

        Parameters:
        - affected_atoms (dict): Dictionary where keys are step numbers and values are lists of atom indices
                                 that were affected in those steps.
        - max_step (int): The maximum step number, used to calculate color gradients for each step.
        """
        try:
            # Iterate over the affected atoms for each step
            for step, atom_indices in affected_atoms.items():
                # Determine the color based on the step number and the maximum step
                color = self.get_color_for_step(step, max_step)

                if color:
                    # Create a unique color name for each step
                    color_name = f"step_{step}_color"

                    # If this color has not been defined in PyMOL yet, define it
                    if color_name not in self.colors_set:
                        self._pymol.cmd.set_color(color_name, color)  # Define color in PyMOL
                        self.colors_set[color_name] = color  # Mark this color as defined

                    # Apply the color to each affected atom
                    for i in atom_indices:
                        if i not in self.colored_atoms:  # Avoid re-coloring the same atom
                            selection = f"index {i}"  # Select the atom by its index
                            self._pymol.cmd.color(color_name, selection)  # Apply the color in PyMOL
                            self.colored_atoms[i] = step  # Track which step colored this atom

                            # Update the PyMOL visualization
                            self._pymolProcess()

                            # Pause briefly to ensure smooth rendering in PyMOL
                            time.sleep(0.01)

            # Calculate and save the percentage of affected atoms
            self.effected_atom_count = round(len(set(self.colored_atoms.keys())) / self.total_atom_count * 100, 1)

            # Save the affected atom percentage to the keeper file
            with open(self.effected_atom_keeper, 'w') as file:
                file.write(f"{self.effected_atom_count}\n")

        except Exception as Err:
            # Catch and print any error that occurs during the process
            print("ERROR:", Err)

    def set_reference_file(self, ref_file_path, effected_atom_keeper):
        """
        Sets the reference file and affected atom data that will be used for comparison
        in future calculations.

        Parameters:
        - ref_file_path (str): Path to the reference data file.
        - effected_atom_keeper (object): Object or data structure that stores information
                                         about affected atoms.
        """
        self.ref_file_path = ref_file_path
        self.effected_atom_keeper = effected_atom_keeper

    def reset_and_update_colors(self, ref_file_path, pert_file_path, threshold=0.05):
        """
        Resets the current color assignments for atoms and recalculates new color mappings
        based on newly provided simulation data. It compares velocities from the reference file
        with perturbed velocities and assigns new colors to affected atoms based on the differences.

        Parameters:
        - ref_file_path (str): Path to the reference velocities or data file.
        - pert_file_path (str): Path to the perturbed velocities or data file.
        - threshold (float): Threshold value for determining affected atoms based on velocity differences.
                             Default is 0.05.
        """
        # Clear the previous data
        self.cumulative_affected_atoms.clear()  # Clear the accumulated affected atoms
        self.colored_atoms.clear()  # Clear previously assigned colors
        self.colors_set.clear()  # Clear defined colors in PyMOL
        self.total_steps = 0  # Reset the step counter

        # Re-assign the reference file and process new velocities
        self.ref_file_path = ref_file_path
        ref_velocities = self.parse_velocities(ref_file_path)
        pert_velocities = self.parse_velocities(pert_file_path)

        # Compare reference and perturbed velocities
        differences = self.compare_velocities(ref_velocities, pert_velocities)

        # Analyze differences using the threshold to determine affected atoms
        affected_atoms = self.find_affected_atoms(differences, threshold=threshold)

        # Start the color update process based on affected atoms
        self.update_colors(affected_atoms, max_step=len(ref_velocities))

        # Set the total atom count if it hasn't been initialized yet
        if self.total_atom_count is None:
            self.total_atom_count = len(ref_velocities)

        # Redraw the visualization in PyMOL
        self._pymolProcess()

"""
if __name__ == "__main__":
    ref_file_path = 'C:\\Users\\law5_\\Desktop\\MDPerTool_GUI\\mdpertool\\output\\1aki_fixed_ph7.4_2024-09-07_20-25-54\\ref_protein_velocities.xml'
    pert_file_path = 'C:\\Users\\law5_\\Desktop\\MDPerTool_GUI\\mdpertool\\output\\1aki_fixed_ph7.4_2024-09-07_20-25-54\\dis_protein_velocities.xml'

    app = QtWidgets.QApplication()
    window = PymolQtWidget()
    window.set_reference_file(ref_file_path=ref_file_path, effected_atom_keeper="out_keep.txt")

    window.loadMolFile("C:\\Users\\law5_\\Desktop\\MDPerTool_GUI\\mdpertool\\output\\1aki_fixed_ph7.4.pdb")
    window.show_protein_in_cartoon_and_atoms()
    window.show()
    window.monitor_live_file(pert_file_path, response_threshold=0.01)
    window.reset_and_update_colors(ref_file_path, pert_file_path)
    #window.update_atom_colors(affected_atoms)
    # Call this function when you want to stop monitoring
    # time.sleep(15)
    # window.stop_monitoring()

    app.exec_()
    window = PymolQtWidget()
    window.show()
"""
