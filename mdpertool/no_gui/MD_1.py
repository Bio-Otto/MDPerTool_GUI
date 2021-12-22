##########################################################################
# IMPORTS
from openmm import app
import openmm as mm
from openmm.unit import femtosecond, picosecond, nanometer, angstrom
from sys import stdout
from openmm import *
from mdtraj.reporters import XTCReporter


class Dissipation_MD_Engine:

    def __init__(self, pdb_path, state_file='state.xml', protein_ff=None, water_ff=None, time_step=None,
                 nonbondedCutoff=None, use_switching_distance=False, switching_distance=None, Device_Index=False,
                 Device_Index_Number=None, dissipation_total_Steps=None, platform_name=None, properties=None,
                 precision=None, CPU_Threads=None, report_interval=1, write_to_dcd=True, dcd_write_period=1,
                 write_to_xtc=False, xtc_write_period=1, dissipated_traj_name='dissipated_traj', output_directory=None):

        print("Simulation parameters preparing for the start ...")
        self.use_switching_distance = use_switching_distance
        self.Device_Index = Device_Index
        self.Device_Index_Number = Device_Index_Number

        if properties is None:
            # --> PROPERTIES SETTINGS FOR OPENCL PLATFORM
            if platform_name == 'OpenCL' and self.Device_Index is True and self.Device_Index_Number is not None:
                if len(self.Device_Index_Number) == 1:
                    properties = {'OpenCLPrecision': '%s' % precision,
                                  'OpenCLDeviceIndex': '%s' % self.Device_Index_Number[0]}
                    print(properties)

                elif len(self.Device_Index_Number) == 2:
                    properties = {'OpenCLPrecision': '%s' % precision,
                                  'OpenCLDeviceIndex': '%s,%s' % (self.Device_Index_Number[0],
                                                                  self.Device_Index_Number[1])}
                    print(properties)
                else:
                    print("Please Control your specified argument for #properties# in OpenCL device")

            if platform_name == 'OpenCL' and self.Device_Index is False:
                properties = {'OpenCLPrecision': '%s' % precision}
                print(properties)

            # --> PROPERTIES SETTINGS FOR CUDA PLATFORM
            if platform_name == 'CUDA' and self.Device_Index == True and self.Device_Index_Number is not None:
                if len(self.Device_Index_Number) == 1:
                    properties = {'CudaPrecision': '%s' % precision, 'DeterministicForces': 'true',
                                  'CudaDeviceIndex': '%s' % self.Device_Index_Number[0]}

                    print(properties)

                elif len(self.Device_Index_Number) == 2:
                    properties = {'CudaPrecision': '%s' % precision, 'DeterministicForces': 'true',
                                  'CudaDeviceIndex': '%s,%s' % (self.Device_Index_Number[0],
                                                                self.Device_Index_Number[1])}

                    print(properties)
                else:
                    print("Please Control your specified argument for #properties# in Cuda device")

            if platform_name == 'CUDA' and self.Device_Index is False:
                properties = {'CudaPrecision': '%s' % precision, 'DeterministicForces': 'true'}
                print(properties)

            # --> PROPERTIES SETTINGS FOR CPU PLATFORM
            if platform_name == 'CPU':
                print("The CPU platform always uses 'mixed' precision.")
                print("Simulation process will use %s Thread(s)" % CPU_Threads)
                properties = {'CpuThreads': '%s' % CPU_Threads}
                precision = 'mixed'

            # --> PROPERTIES SETTINGS FOR REFERENCE PLATFORM
            if platform_name == 'Reference':
                print("The Reference platform always uses 'double' precision.")
                properties = None
                precision = 'double'


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

        self.state_file = state_file
        self.switching_distance = switching_distance * angstrom
        self.time_step = time_step * femtosecond
        self.nonbondedCutoff = nonbondedCutoff * angstrom
        self.dissipation_total_Steps = dissipation_total_Steps
        self.platform_name = platform_name
        self.properties = properties
        self.report_interval = report_interval
        self.write_to_dcd = write_to_dcd
        self.dcd_write_period = dcd_write_period
        self.write_to_xtc = write_to_xtc
        self.xtc_write_period = xtc_write_period
        self.output_directory = output_directory

        # we'll just take the topology from here...
        pdb = app.PDBFile(self.pdb_path)
        topology = pdb.topology

        print('Forcefield parameters loading to the simulation system ...')
        # forcefield = app.ForceField('amber96.xml', 'tip3p.xml')
        forcefield = app.ForceField(self.protein_ff, self.water_ff)

        print('Constructing an OpenMM System')
        # self.system = forcefield.createSystem(topology, nonbondedMethod=app.PME, nonbondedCutoff=self.nonbondedCutoff,
        #                                       constraints=None, rigidWater=True, ewaldErrorTolerance=1e-5)
        self.system = forcefield.createSystem(topology, nonbondedMethod=app.PME, nonbondedCutoff=self.nonbondedCutoff,
                                              constraints=None, rigidWater=True, ewaldErrorTolerance=1e-5)
        # print("USES OR NOT: %s" %self.system.usesPeriodicBoundaryConditions())

        if self.use_switching_distance:
            print("System will use Switching Distance")
            nonbonded = [f for f in self.system.getForces() if isinstance(f, NonbondedForce)][0]
            nonbonded.setUseSwitchingFunction(use=True)
            nonbonded.setSwitchingDistance(self.switching_distance)
            # nonbonded.setUseDispersionCorrection(False)
            # print("Dispertion Correction Not working")

            # nonbonded.setReciprocalSpaceForceGroup(1)

        # state = mm.XmlSerializer.deserialize(open(self.state_file).read())

        integrator = VerletIntegrator(self.time_step)
        integrator.setConstraintTolerance(1e-8)

        # let's specify our simulation platform again
        platform = mm.Platform.getPlatformByName(self.platform_name)

        # ok now let's do some simulation using this restraint
        simulation = app.Simulation(topology, self.system, integrator, platform, self.properties)

        simulation.loadState(self.state_file)
        final_state = simulation.context.getState(getVelocities=True, getPositions=True)
        positions = final_state.getPositions()
        velocities = final_state.getVelocities()
        simulation.context.setPositions(positions)
        simulation.context.computeVirtualSites()
        simulation.context.setVelocities(velocities)
        print("USES OR NOT: %s" % self.system.usesPeriodicBoundaryConditions())
        # simulation.context.setTime(0)

        if self.report_interval > self.dissipation_total_Steps:
            self.report_interval = int(self.dissipation_total_Steps / 2)
            print("The number of report steps has been adjusted to %s, because the number of report steps exceeds "
                  "the total number of steps." % self.report_interval)

        if self.write_to_dcd:
            print("Saving DCD File for every %s period" % self.dcd_write_period)
            DCD_file_path = os.path.join(self.output_directory, '%s.dcd' % dissipated_traj_name)
            simulation.reporters.append(app.DCDReporter(DCD_file_path, self.dcd_write_period))

        if self.write_to_xtc:
            print("Saving XTC File for every %s period" % self.xtc_write_period)
            XTC_file_path = os.path.join(self.output_directory, '%s.xtc' % dissipated_traj_name)
            simulation.reporters.append(XTCReporter(XTC_file_path, self.xtc_write_period))

        simulation.reporters.append(
            app.StateDataReporter(stdout, self.report_interval, step=True, time=True, potentialEnergy=True,
                                  kineticEnergy=True, totalEnergy=True, temperature=True, progress=True, volume=True,
                                  density=True, remainingTime=True, speed=True, totalSteps=self.dissipation_total_Steps,
                                  separator='\t'))

        # # set up reporters so we can see what's going on...
        # simulation.reporters.append(app.PDBReporter(self.multi_pdb_name, self.snapshot_period))

        simulation.step(self.dissipation_total_Steps)
