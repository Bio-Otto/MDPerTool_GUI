from simtk.openmm.app import StateDataReporter
from io import StringIO
import time
import queue
import threading
import itertools
import tokenize
import os
import pystache
import numpy as np
# from chaco.api import Plot, ArrayPlotData, PlotAxis, VPlotContainer

import pyqtgraph as pg
#
#
# ##############################################################################
# # Functions
# ##############################################################################
#
#
# def queue_reporter_factory(queue):
#     """Factory function that returns a dynamically defined OpenMM
#     reporter class which reports by sending dicts down a synchronous queue"""
#
#     class QueueStateDataReporter(StateDataReporter):
#         """Subclass of StateDataReporter sends its results down a synchronous
#         Queue, as opposed to printing them to a file-like object
#         """
#
#         def __init__(self, file, *args, **kwargs):
#             with open(os.devnull, 'w') as f:
#                 # send in a fake file
#                 super(QueueStateDataReporter, self).__init__(f, *args, **kwargs)
#
#             # this is where we'll store the names of the fields that
#             # are being reported in
#             self._headers = []
#
#         def report(self, simulation, state):
#             was_initialized = self._hasInitialized
#
#             # spoof the file-like object with a string buffer
#             self._out = StringIO()
#             super(QueueStateDataReporter, self).report(simulation, state)
#
#             if not was_initialized:
#                 # the first report has two lines on it -- we want to
#                 # look at the first, as it contains the headers
#                 initial, line = self._out.getvalue().split(os.linesep, 1)
#                 headers = initial.strip().split(self._separator)
#                 # filter out some extra quotation marks and comment
#                 # characters
#                 self._headers = [e.strip('#"\'') for e in headers]
#             else:
#                 line = self._out.getvalue()
#
#             # split the line based on whatever separator we know
#             # that the parent was using, and then cast to float
#             values = map(float, line.strip().split(self._separator))
#             msg = dict(zip(self._headers, values))
#             queue.put(msg)
#
#     return QueueStateDataReporter
#
#
# def run_openmm_script(code, queue):
#     """Run an OpenMM script, dynamiclly redefining a StateDataReporter
#     to capture the output, sending it into this queue
#     """
#
#     def fix_code():
#         """Replae the token 'StateDataReporter' with
#         '__queue_reporter_factory(__queue)'
#         Also, we make sure that the sentenel signal (None) is sent
#         down the queue at the very end of the script
#         """
#         itoks = tokenize.generate_tokens(StringIO(code).readline)
#
#         def run():
#             for toktype, toktext, (srow, scol), (erow, ecol), line in itoks:
#                 if toktext == 'StateDataReporter':
#                     toktext = '__queue_reporter_factory(__queue)'
#                 yield (toktype, toktext, (srow, scol), (erow, ecol), line)
#
#         return tokenize.untokenize(run()) + '__queue.put(None)'
#
#     try:
#         code = fix_code()
#     except tokenize.TokenError:
#         raise ValueError('The script has a syntax error!')
#
#     exec(code, {'__queue': queue, '__queue_reporter_factory': queue_reporter_factory})
#
#
# def queue_consumer(q):
#     """Main loop for a thread that consumes the messages from the queue
#     and plots them"""
#
#     status = 'Running...'
#
#     while True:
#         try:
#             msg = q.get_nowait()
#             if msg is None:
#                 print("None")
#                 break
#             update_plot(msg)
#         except queue.Empty:
#             time.sleep(0.1)
#
#     status = 'Done'
# def update_plot(msg):
#     """Add data points from the message to the plots
#     Paramters
#     ---------
#     msg : dict
#         This is the message sent over the Queue from the script
#     """
#     print(msg)
#     # if not plots_created:
#     #     self.create_plots(msg.keys())
#     #     self.plots_created = True
#     #
#     # for k, v in msg.iteritems():
#     #     current = self.plotdata.get_data(k)
#     #     self.plotdata.set_data(k, np.r_[current, v])
#
renderer = pystache.Renderer()
template = pystache.parse(u'''

from simtk.openmm import app
from simtk.openmm.app import *
from simtk.openmm.app import PME, NoCutoff, Ewald, CutoffPeriodic, CutoffNonPeriodic, HBonds, HAngles, AllBonds
import simtk.openmm as mm
from simtk.unit import femtosecond, picosecond, nanometer, kelvin, angstrom, atmospheres
from sys import stdout
from apply_pdbfixer import fix_pdb
from simtk.openmm import *
from mdtraj.reporters import XTCReporter

print('pdb file fixing and preparing for simulation ...')
fixed_pdb_name = fix_pdb('C:/Users/HIbrahim/Desktop/MDPERTOOL_v01/Download/1gg1.pdb')

pdb = app.PDBFile(fixed_pdb_name)

box = pdb.topology.getUnitCellDimensions()

print('Modeller of pdb file is preparing ...')
modeller = app.Modeller(pdb.topology, pdb.positions)
modeller.topology.setUnitCellDimensions(box)

print('Forcefield parameters loading to the simulation system ...')
forcefield = app.ForceField('amber96.xml', 'tip3p.xml')

print('Adding missing hydrogens to the model ...')
modeller.addHydrogens(forcefield)

print('Adding solvent (both water and ions) to the model to fill a rectangular box ...')
modeller.addSolvent(forcefield, model='tip3p', padding=15*angstrom)

print('Constructing an OpenMM System')
system = forcefield.createSystem(modeller.topology, nonbondedMethod=PME, nonbondedCutoff=12.0*angstrom, constraints=None, rigidWater=True, ewaldErrorTolerance=0.005)


integrator = mm.LangevinIntegrator(310*unit.kelvin, 1.0/unit.picoseconds, 2.0*unit.femtoseconds)
integrator.setConstraintTolerance(0.00001)

platform = mm.Platform.getPlatformByName('OpenCL')
properties = {'OpenCLPrecision': 'single'}
simulation = app.Simulation(modeller.topology, system, integrator, platform, properties)
simulation.context.setPositions(modeller.positions)

print('Minimizing...')
simulation.minimizeEnergy()

simulation.context.setVelocitiesToTemperature(300*unit.kelvin)
print('Equilibrating...')
simulation.step(100)

simulation.reporters.append(app.DCDReporter('trajectory.dcd', 1000))
simulation.reporters.append(StateDataReporter(stdout, 10, step=True, potentialEnergy=True, temperature=True, progress=True, remainingTime=True, speed=True, totalSteps=1000))

print('Running Production...')
simulation.step(1000)
print('Done!')
''')

contents = renderer.render(template)


#
# q = queue.Queue()
# # start up two threads. the first, t1, will run the script
# # and place the statedata into the queue
# # the second will remove elements from the queue and update the
# # plots in the UI
# t1 = threading.Thread(target=run_openmm_script,
#                       args=(contents, q))
# t2 = threading.Thread(target=queue_consumer, args=(q,))
# t1.start()
# t2.start()


##############################################################################
# Functions
##############################################################################


def queue_reporter_factory(queue):
    """Factory function that returns a dynamically defined OpenMM
    reporter class which reports by sending dicts down a synchronous queue"""

    class QueueStateDataReporter(StateDataReporter):
        """Subclass of StateDataReporter sends its results down a synchronous
        Queue, as opposed to printing them to a file-like object
        """

        def __init__(self, file, *args, **kwargs):
            with open(os.devnull, 'w') as f:
                # send in a fake file
                super(QueueStateDataReporter, self).__init__(f, *args, **kwargs)

            # this is where we'll store the names of the fields that
            # are being reported in
            self._headers = []

        def report(self, simulation, state):
            was_initialized = self._hasInitialized

            # spoof the file-like object with a string buffer
            self._out = StringIO()
            super(QueueStateDataReporter, self).report(simulation, state)

            if not was_initialized:
                # the first report has two lines on it -- we want to
                # look at the first, as it contains the headers
                print(self._out.getvalue())
                initial, line = self._out.getvalue().split('\n', 1)
                headers = initial.strip().split(',')
                # filter out some extra quotation marks and comment
                # characters
                self._headers = [e.strip('#"\'') for e in headers]
            else:
                line = self._out.getvalue()

            t = [e.strip('%"\'') for e in line.strip().split(',')]

            # split the line based on whatever separator we know
            # that the parent was using, and then cast to float
            # values = map(float, line.strip().split(self._separator))
            msg = dict(zip(self._headers, t))
            queue.put(msg)

    return QueueStateDataReporter


def run_openmm_script(code, queue):
    """Run an OpenMM script, dynamiclly redefining a StateDataReporter
    to capture the output, sending it into this queue
    """

    def fix_code():
        """Replae the token 'StateDataReporter' with
        '__queue_reporter_factory(__queue)'
        Also, we make sure that the sentenel signal (None) is sent
        down the queue at the very end of the script
        """
        itoks = tokenize.generate_tokens(StringIO(code).readline)

        def run():
            for toktype, toktext, (srow, scol), (erow, ecol), line in itoks:
                if toktext == 'StateDataReporter':
                    toktext = '__queue_reporter_factory(__queue)'
                yield toktype, toktext, (srow, scol), (erow, ecol), line

        return tokenize.untokenize(run()) + '__queue.put(None)'

    try:
        code = fix_code()
    except tokenize.TokenError:
        raise ValueError('The script has a syntax error!')

    exec(code, {'__queue': queue, '__queue_reporter_factory': queue_reporter_factory})


def chaco_scatter(dataview, x_name, y_name, x_label=None, y_label=None,
                  color=None):
    """Utility function to build a chaco scatter plot
    """
    plot = dataview
    print("PLOT:%s" %plot)
    # plot = Plot(dataview)
    # plot.plot((x_name, y_name), type="scatter", marker='dot', color=color)
    #
    # if x_label is None:
    #     x_label = x_name
    # if y_label is None:
    #     y_label = y_name
    # x_axis = PlotAxis(mapper=plot.x_mapper, orientation='bottom', title=x_label)
    # y_axis = PlotAxis(mapper=plot.y_mapper, orientation='left', title=y_label)
    # plot.underlays.append(x_axis)
    # plot.underlays.append(y_axis)
    return plot


##############################################################################
# Classes
##############################################################################


class OpenMMScriptRunner:
    plots = Container
    plots_created = bool
    openmm_script_code = str
    status = str

    def __init__(self, script):
        super(OpenMMScriptRunner, self).__init__()

        self._plots_created = False
        self.openmm_script_code = script
        q = queue.Queue()
        # start up two threads. the first, t1, will run the script and place the statedata into the queue
        # the second will remove elements from the queue and update the plots in the UI
        t1 = threading.Thread(target=run_openmm_script, args=(self.openmm_script_code, q))
        t2 = threading.Thread(target=self.queue_consumer, args=(q,))
        t1.start()
        t2.start()

    def queue_consumer(self, q):
        """Main loop for a thread that consumes the messages from the queue
        and plots them"""

        self.status = 'Running...'
        while True:

            try:
                msg = q.get_nowait()

                if msg is None:
                    print("geldiii None")
                    break
                self.update_plot(msg)
            except queue.Empty:
                time.sleep(0.1)

        self.status = 'Done'

    def create_plots(self, keys):
        """Create the plots
        Paramters
        ---------
        keys : list of strings
            A list of all of the keys in the msg dict. This should be something
            like ['Step', 'Temperature', 'Potential Energy']. We'll create the
            ArrayPlotData container in which each of these timeseries will
            get put.
        """
        print("KEYS: %s" % keys)
        # self.plots = Container
        # # this looks cryptic, but it is equivalent to
        # # ArrayPlotData(a=[], b=[], c=[])
        # # if the keys are a,b,c. This just does it for all of the keys.
        # self.plotdata = (**dict(zip(keys, [[]] * len(keys))))
        #
        # # figure out which key will be the x axis
        # if 'Step' in keys:
        #     x = 'Step'
        # elif 'Time (ps)' in keys:
        #     x = 'Time (ps)'
        # else:
        #     raise ValueError('The reporter published neither the step nor time'
        #                      'count, so I don\'t know what to plot on the x-axis!')
        #
        # colors = itertools.cycle(['blue', 'green', 'silver', 'pink', 'lightblue',
        #                           'red', 'darkgray', 'lightgreen', ])
        # for y in filter(lambda y: y != x, keys):
        #     self.plots.add(chaco_scatter(self.plotdata, x_name=x, y_name=y,
        #                                  color=colors.next()))

    def update_plot(self, msg):
        """Add data points from the message to the plots
        Paramters
        ---------
        msg : dict
            This is the message sent over the Queue from the script
        """
        print("MESSAGE: %s" % msg)
        # if not self.plots_created:
        #     self.create_plots(msg.keys())
        #     self.plots_created = True
        #
        # for k, v in msg.iteritems():
        #     current = self.plotdata.get_data(k)
        #     self.plotdata.set_data(k, np.r_[current, v])


runner = OpenMMScriptRunner(contents)
# runner.configure_traits()