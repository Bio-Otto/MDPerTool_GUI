##########################################################################
# IMPORTS
from simtk.openmm import app
import simtk.openmm as mm
from simtk.unit import femtosecond, picosecond, nanometer, kelvin, angstrom, atmospheres
from sys import stdout
from .apply_pdbfixer import fix_pdb
from simtk.openmm import *
from mdtraj.reporters import XTCReporter


class Classic_MD_Engine:
    def __init__(self, pdb_path, protein_ff='amber96', water_ff='tip3p', time_step=2.0, nonbondedCutoff=12.0,
                 switching_distance = 10.0, water_padding=15, Device_Index=False, Device_Index_Number=1,
                 total_Steps=300000, temp=310, platform_name='OpenCL', properties=None, precision='single',
                 friction_cofficient=1.0, minimize=True, minimize_steps=5000, CPU_Threads=2, equilibrate=True,
                 equilibration_step=500, report_interval=500, write_system_xml=False, system_file_name='system.xml',
                 state_file_name='state.xml', last_pdb_filename='last.pdb', write_to_dcd=False, dcd_write_period=50,
                 write_to_xtc=False, xtc_write_period=50):

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


        self.time_step = time_step * femtosecond
        self.nonbondedCutoff = nonbondedCutoff * angstrom
        self.switching_distance = switching_distance * angstrom
        self.water_padding = water_padding * angstrom
        self.total_Steps = total_Steps
        self.temp = temp * kelvin
        self.platform_name = platform_name
        self.properties = properties
        self.friction_cofficient = friction_cofficient / picosecond
        self.minimize = minimize
        self.minimize_steps = minimize_steps
        self.equilibrate = equilibrate
        self.equilibration_step = equilibration_step
        self.report_interval = report_interval
        self.write_system_xml = write_system_xml
        self.system_file_name = system_file_name
        self.state_file_name = state_file_name
        self.last_pdb_filename = last_pdb_filename

        self.write_to_dcd = write_to_dcd
        self.dcd_write_period = dcd_write_period
        self.write_to_xtc = write_to_xtc
        self.xtc_write_period = xtc_write_period


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

        ## FIXING PDB
        print('pdb file fixing and preparing for simulation ...')
        fixed_pdb_name = fix_pdb(self.pdb_path)

        print('Loading pdb to simulation engine ...')
        pdb = app.PDBFile(fixed_pdb_name)

        box = pdb.topology.getUnitCellDimensions()

        print('Modeller of pdb file is preparing ...')
        modeller = mm.app.Modeller(pdb.topology, pdb.positions)
        modeller.topology.setUnitCellDimensions(box)

        print('Forcefield parameters loading to the simulation system ...')
        forcefield = app.ForceField(self.protein_ff, self.water_ff)

        print('Adding missing hydrogens to the model ...')
        modeller.addHydrogens(forcefield)

        print('Adding solvent (both water and ions) to the model to fill a rectangular box ...')
        modeller.addSolvent(forcefield, model=self.water_model, padding=self.water_padding)

        print('Constructing an OpenMM System')
        self.system = forcefield.createSystem(modeller.topology, nonbondedMethod=app.PME,
                                              nonbondedCutoff=self.nonbondedCutoff, constraints=None,
                                              rigidWater=True,
                                              ewaldErrorTolerance=0.005)

        # self.system.addForce(mm.AndersenThermostat(310 * unit.kelvin, self.friction_cofficient))

        self.system.addForce(mm.MonteCarloBarostat(1 * atmospheres, self.temp, 25))

        nonbonded = [f for f in self.system.getForces() if isinstance(f, NonbondedForce)][0]
        nonbonded.setUseSwitchingFunction(use=True)
        nonbonded.setSwitchingDistance(self.switching_distance)
        nonbonded.setUseDispersionCorrection(True)

        print('Creating a LangevinIntegrator.')
        integrator = mm.LangevinIntegrator(self.temp, self.friction_cofficient, self.time_step)

        platform = mm.Platform.getPlatformByName(self.platform_name)

        simulation = app.Simulation(modeller.topology, self.system, integrator, platform, self.properties)
        simulation.context.setPositions(modeller.positions)

        simulation.context.computeVirtualSites()

        if minimize:
            print('Minimizing...')
            simulation.minimizeEnergy(maxIterations=self.minimize_steps)
            print("Minimization done, the energy is", simulation.context.getState(getEnergy=True).getPotentialEnergy())
            positions = simulation.context.getState(getPositions=True).getPositions()
            print("Minimized geometry is written to 'minimized.pdb'")
            app.PDBFile.writeModel(modeller.topology, positions, open('minimized.pdb', 'w'), keepIds=True)

        simulation.context.setVelocitiesToTemperature(self.temp)

        if self.equilibrate:
            print('Equilibrating...')
            simulation.step(self.equilibration_step)

        if self.report_interval > self.total_Steps:
            self.report_interval = int(self.total_Steps / 2)
            print("The number of report steps has been adjusted to %s, because the number of report steps exceeds "
                  "the total number of steps." % self.report_interval)

        if self.write_to_dcd:
            print("Saving DCD File for every %s period" % self.dcd_write_period)
            simulation.reporters.append(app.DCDReporter('trajectory.dcd', self.dcd_write_period))

        if self.write_to_xtc:
            print("Saving XTC File for every %s period" % self.xtc_write_period)
            simulation.reporters.append(XTCReporter('trajectory.xtc', self.xtc_write_period))

        simulation.reporters.append(app.StateDataReporter(stdout, self.report_interval, step=True,
                                                          time=True, potentialEnergy=True, kineticEnergy=True,
                                                          totalEnergy=True, temperature=True, progress=True,
                                                          remainingTime=True, speed=True, volume=True, density=True,
                                                          totalSteps=self.total_Steps, separator='\t'))

        if self.write_system_xml:
            print("Serializing the system")
            serial = mm.XmlSerializer.serializeSystem(self.system)
            with open('equilubrate_system.xml', 'w') as f:
                f.write(serial)

            # Print out information about the platform if we wan

        print('Running Production...')
        simulation.step(self.total_Steps)  # 300000
        print('Done!')

        lastpositions = simulation.context.getState(getPositions=True).getPositions()

        self.last_pdb = app.PDBFile.writeFile(modeller.topology, lastpositions, open(self.last_pdb_filename, 'w'),
                                              keepIds=True)

        print("###################################")
        # Now we will save a serialization of this simulation into OpenMM's native XML format
        # We can re-initialize the system later for further simulations without all of the bothersome set-up by
        # loading these files!

        self.state = simulation.context.getState(getPositions=True, getVelocities=True)

        # system.xml contains all of the force field parameters

        with open(self.system_file_name, 'w') as f:
            system_xml = mm.XmlSerializer.serialize(self.system)
            f.write(system_xml)
        # integrator.xml contains the confiruation for the integrator, RNG seed
        with open('integrator.xml', 'w') as f:
            integrator_xml = mm.XmlSerializer.serialize(integrator)
            f.write(integrator_xml)
            # state.xml contains positions, velocities, forces, the barostat
        with open(self.state_file_name, 'w') as f:
            f.write(mm.XmlSerializer.serialize(self.state))

        simulation.context.setTime(0)

        

# Reference_MD_Engine('2j0w.pdb', minimize=True, minimize_steps=50, platform_name='OpenCL', water_padding=10)