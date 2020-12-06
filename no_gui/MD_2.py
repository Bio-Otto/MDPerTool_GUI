##########################################################################
# IMPORTS
from simtk.openmm import app
import simtk.openmm as mm
from simtk.unit import femtosecond, picosecond, nanometer, angstrom
from sys import stdout
from simtk.openmm import *
from mdtraj.reporters import XTCReporter


class Reference_MD_Engine:

    def __init__(self, pdb_path, state_file='state.xml', protein_ff='amber96', water_ff='tip3p', time_step=1.0,
                 nonbondedCutoff=12.0, switching_distance= 10.0, Device_Index=False, Device_Index_Number=1,
                 reference_total_Steps=1000, platform_name='OpenCL', properties=None, precision='single',
                 CPU_Threads=2, report_interval=1, write_to_dcd=True, dcd_write_period=1, write_to_xtc=False,
                 xtc_write_period=1, undissipated_traj_name='undissipated_traj'):

        print("Simulation parameters preparing for the start ...")
        self.Device_Index = Device_Index
        self.Device_Index_Number = Device_Index_Number

        if properties is None:
            if platform_name == 'OpenCL' and self.Device_Index == True:
                properties = {'OpenCLPrecision': '%s' % precision, 'OpenCLDeviceIndex': '%s' % self.Device_Index_Number}

            if platform_name == 'OpenCL' and self.Device_Index == False:
                properties = {'OpenCLPrecision': '%s' % precision}

            if platform_name == 'CUDA' and self.Device_Index == True:
                properties = {'CudaPrecision': '%s' % precision, 'CudaDeviceIndex': '%s' % self.Device_Index_Number}

            if platform_name == 'CPU':
                print("The CPU platform always uses 'mixed' precision.")
                print("Simulation process will use %s Thread(s)" % CPU_Threads)
                properties = {'CpuThreads': '%s' % CPU_Threads}
                precision = 'mixed'

            if platform_name == 'Reference':
                print("The Reference platform always uses 'double' precision.")
                properties = None
                precision = 'double'
        print("System will use '%s' Platform with '%s' Precision" % (platform_name, precision))

        self.pdb_path = pdb_path
        self.protein_ff = protein_ff + '.xml'

        if protein_ff == 'charmm36':
            print("Protein FF: %s" % self.protein_ff)
            self.water_ff = '%s/%s.xml' % (protein_ff, water_ff)
            print("Water FF: %s" % self.water_ff)
        else:
            self.water_ff = '%s.xml' % water_ff
            print("Protein FF: %s" % self.protein_ff)
            print("Water FF: %s" % self.water_ff)

        ## SOLUTION WATER MODEL
        if len((self.water_ff.split('/')[-1]).split('.')[0]) <= 5:
            self.water_model = (self.water_ff.split('/')[-1]).split('.')[0]
            print("Water Model for Solution 1: %s" % self.water_model)

        elif (self.water_ff.split('/')[-1]).split('.')[0] == 'tip4pew':
            self.water_model = (self.water_ff.split('/')[-1]).split('.')[0]
            print("Water Model for Solution 2: %s" % self.water_model)

        else:
            self.water_model = (self.water_ff.split('/')[-1]).split('.')[0][0:5]
            print("Water Model for Solution 3: %s" % self.water_model)

        self.switching_distance = switching_distance * angstrom
        self.state_file = state_file
        self.time_step = time_step * femtosecond
        self.nonbondedCutoff = nonbondedCutoff * angstrom
        self.reference_total_Steps = reference_total_Steps
        self.platform_name = platform_name
        self.properties = properties
        self.report_interval = report_interval
        self.write_to_dcd = write_to_dcd
        self.dcd_write_period = dcd_write_period
        self.write_to_xtc = write_to_xtc
        self.xtc_write_period = xtc_write_period


        # we'll just take the topology from here...
        pdb = app.PDBFile(self.pdb_path)
        topology = pdb.topology

        print('Forcefield parameters loading to the simulation system ...')
        # forcefield = app.ForceField('amber96.xml', 'tip3p.xml')
        forcefield = app.ForceField(self.protein_ff, self.water_ff)

        print('Constructing an OpenMM System')
        self.system = forcefield.createSystem(topology, nonbondedMethod=app.PME, nonbondedCutoff=self.nonbondedCutoff,
                                              constraints=None, rigidWater=True, ewaldErrorTolerance=5e-5)

        nonbonded = [f for f in self.system.getForces() if isinstance(f, NonbondedForce)][0]
        nonbonded.setUseSwitchingFunction(use=True)
        nonbonded.setSwitchingDistance(self.switching_distance)
        # nonbonded.setUseDispersionCorrection(True)
        # nonbonded.setReciprocalSpaceForceGroup(1)

        state = mm.XmlSerializer.deserialize(open(self.state_file).read())

        integrator = VerletIntegrator(self.time_step)

        # let's specify our simulation platform again
        platform = mm.Platform.getPlatformByName(self.platform_name)

        # ok now let's do some simulation using this restraint
        simulation = app.Simulation(topology, self.system, integrator, platform, self.properties)
        simulation.context.setState(state)
        simulation.context.setTime(0)

        simulation.context.computeVirtualSites()

        if self.report_interval > self.reference_total_Steps:
            self.report_interval = int(self.reference_total_Steps / 2)
            print("The number of report steps has been adjusted to %s, because the number of report steps exceeds "
                  "the total number of steps." % self.report_interval)

        if self.write_to_dcd:
            print("Saving DCD File for every %s period" % self.dcd_write_period)
            simulation.reporters.append(app.DCDReporter('%s.dcd' % undissipated_traj_name, self.dcd_write_period))

        if self.write_to_xtc:
            print("Saving XTC File for every %s period" % self.xtc_write_period)
            simulation.reporters.append(XTCReporter('%s.xtc' % undissipated_traj_name, self.xtc_write_period))

        simulation.reporters.append(app.StateDataReporter(stdout, self.report_interval, step=True, time=True, potentialEnergy=True,
                                  kineticEnergy=True, totalEnergy=True, temperature=True, progress=True, volume=True,
                                  density=True, remainingTime=True, speed=True, totalSteps=self.reference_total_Steps,
                                  separator='\t'))

        # # set up reporters so we can see what's going on...
        # simulation.reporters.append(app.PDBReporter(self.multi_pdb_name, self.snapshot_period))

        simulation.step(self.reference_total_Steps)
