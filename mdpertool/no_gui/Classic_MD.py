##########################################################################
# IMPORTS
from openmm import app
import openmm as mm
from openmm.unit import femtosecond, picosecond, nanometer, kelvin, angstrom, atmospheres
from sys import stdout
from .apply_pdbfixer import fix_pdb
from openmm import *
from mdtraj.reporters import XTCReporter
from .optimizepme import *


class Classic_MD_Engine:
    def __init__(self, pdb_path, protein_ff=None, water_ff=None, time_step=2.0, nonbondedCutoff=12.0,
                 use_switching_distance=True,
                 switching_distance=10.0, water_padding=None, Device_Index=False, Device_Index_Number=None,
                 total_Steps=300000, temp=310, platform_name=None, properties=None, precision=None,
                 friction_cofficient=91.0, minimize=True, minimize_steps=5000, just_minimize=False, CPU_Threads=None,
                 equilibrate=True, equilibration_step=500, report_interval=500, write_system_xml=False,
                 system_file_name='system.xml', state_file_name='state.xml', last_pdb_filename='last.pdb',
                 write_to_dcd=False, dcd_write_period=50, write_to_xtc=False, xtc_write_period=50,
                 output_directory=None, user_selected_res=None):

        print("Simulation parameters preparing for the start ...")
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
                    properties = {'CudaPrecision': '%s' % precision,
                                  'CudaDeviceIndex': '%s' % self.Device_Index_Number[0], 'DeterministicForces': 'True'}

                    print(properties)

                elif len(self.Device_Index_Number) == 2:
                    properties = {'CudaPrecision': '%s' % precision,
                                  'CudaDeviceIndex': '%s,%s' % (self.Device_Index_Number[0],
                                                                self.Device_Index_Number[1]),
                                  'DeterministicForces': 'True'}

                    print(properties)
                else:
                    print("Please Control your specified argument for #properties# in Cuda device")

            if platform_name == 'CUDA' and self.Device_Index is False:
                properties = {'CudaPrecision': '%s' % precision, 'DeterministicForces': 'True'}
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
        self.use_switching_distance = use_switching_distance
        self.switching_distance = switching_distance * angstrom
        self.water_padding = water_padding * angstrom
        self.total_Steps = total_Steps
        self.temp = temp * kelvin
        self.platform_name = platform_name
        self.properties = properties
        self.friction_cofficient = friction_cofficient / picosecond
        self.minimize = minimize
        self.minimize_steps = minimize_steps
        self.just_minimize = just_minimize
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
        self.output_directory = output_directory

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
        fixed_pdb_name = fix_pdb(self.pdb_path, self.output_directory)

        print('Loading pdb to simulation engine ...')
        pdb = app.PDBFile(fixed_pdb_name)

        print("Perturbed Residue that you chose is controlling")
        residue_control = [r for r in pdb.topology.residues()]  # build a list of residues

        # print(residue_control)
        pert_res_len = len(user_selected_res)
        pert_len_must = 0
        for index in range(len(residue_control)):
            ind = residue_control[index].id
            nm = residue_control[index].name
            resi = str(nm) + str(ind)
            if resi in user_selected_res:
                print("TOPOLOGY INCLUDING RESIDUE %s" % resi)
                pert_len_must += 1

        if pert_res_len != pert_len_must:
            print("OHH NO! YOUR SELECTED RESIDUES NOT IN TOPOLOGY")
            import shutil
            try:
                shutil.rmtree(output_directory)
            except OSError as e:
                print("Error: %s : %s" % (output_directory, e.strerror))

            sys.exit(0)

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
                                              rigidWater=None,
                                              ewaldErrorTolerance=0.00001)

        # self.system.addForce(mm.AndersenThermostat(310 * unit.kelvin, self.friction_cofficient))

        self.system.addForce(mm.MonteCarloBarostat(1.0 * atmospheres, self.temp, 25))

        if self.use_switching_distance:
            print("System will use Switching Distance")
            nonbonded = [f for f in self.system.getForces() if isinstance(f, NonbondedForce)][0]
            nonbonded.setUseSwitchingFunction(use=True)
            nonbonded.setSwitchingDistance(self.switching_distance)
            nonbonded.setUseDispersionCorrection(False)

        print('Creating a LangevinIntegrator.')
        integrator = mm.LangevinMiddleIntegrator(self.temp, self.friction_cofficient, self.time_step)

        platform = mm.Platform.getPlatformByName(self.platform_name)

        simulation = app.Simulation(modeller.topology, self.system, integrator, platform, self.properties)
        simulation.context.setPositions(modeller.positions)

        simulation.context.computeVirtualSites()

        if minimize:
            print('Minimizing...')
            minimize_pdb_path = os.path.join(self.output_directory, 'minimized.pdb')
            simulation.minimizeEnergy(maxIterations=self.minimize_steps)
            print("Minimization done, the energy is", simulation.context.getState(getEnergy=True).getPotentialEnergy())
            positions = simulation.context.getState(getPositions=True).getPositions()
            print("Minimized geometry is written to 'minimized.pdb'")

            app.PDBFile.writeModel(modeller.topology, positions, open(minimize_pdb_path, 'w'), keepIds=True)

        simulation.context.setVelocitiesToTemperature(self.temp)

        if not self.just_minimize:
            if self.equilibrate:
                print('Equilibrating...')
                simulation.step(self.equilibration_step)

            simulation.context.setTime(0)

            if self.report_interval > self.total_Steps:
                self.report_interval = int(self.total_Steps / 2)
                print("The number of report steps has been adjusted to %s, because the number of report steps exceeds "
                      "the total number of steps." % self.report_interval)

            if self.write_to_dcd:
                DCD_file_path = os.path.join(self.output_directory, 'trajectory.dcd')
                print("Saving DCD File for every %s period" % self.dcd_write_period)
                simulation.reporters.append(app.DCDReporter(DCD_file_path, self.dcd_write_period))

            if self.write_to_xtc:
                XTC_file_path = os.path.join(self.output_directory, 'trajectory.xtc')
                print("Saving XTC File for every %s period" % self.xtc_write_period)
                simulation.reporters.append(XTCReporter(XTC_file_path, self.xtc_write_period))

            simulation.reporters.append(app.StateDataReporter(stdout, self.report_interval, step=True,
                                                              time=True, potentialEnergy=True, kineticEnergy=True,
                                                              totalEnergy=True, temperature=True, progress=True,
                                                              remainingTime=True, speed=True, volume=True, density=True,
                                                              totalSteps=self.total_Steps, separator='\t'))

            if self.write_system_xml:
                System_XML_file_path = os.path.join(self.output_directory, 'equilubrate_system.xml')
                print("Serializing the system")
                serial = mm.XmlSerializer.serializeSystem(self.system)
                with open(System_XML_file_path, 'w') as f:
                    f.write(serial)

                # Print out information about the platform if we wan

            print('Running Production...')
            simulation.step(self.total_Steps)  # 300000
            print('Done!')

            lastpositions = simulation.context.getState(getPositions=True).getPositions()
            last_pdb_file_path = os.path.join(self.output_directory, self.last_pdb_filename)
            self.last_pdb = app.PDBFile.writeFile(modeller.topology, lastpositions, open(last_pdb_file_path, 'w'),
                                                  keepIds=True)

            print("###################################")
            # Now we will save a serialization of this simulation into OpenMM's native XML format
            # We can re-initialize the system later for further simulations without all of the bothersome set-up by
            # loading these files!

            self.state = simulation.context.getState(getPositions=True, getVelocities=True)

            # system.xml contains all of the force field parameters
            self.system_file_path = os.path.join(self.output_directory, self.system_file_name)
            with open(self.system_file_path, 'w') as f:
                system_xml = mm.XmlSerializer.serialize(self.system)
                f.write(system_xml)

            # integrator.xml contains the confiruation for the integrator, RNG seed
            self.integrator_file_path = os.path.join(self.output_directory, 'integrator.xml')
            with open(self.integrator_file_path, 'w') as f:
                integrator_xml = mm.XmlSerializer.serialize(integrator)
                f.write(integrator_xml)

            # state.xml contains positions, velocities, forces, the barostat
            self.state_file_path = os.path.join(self.output_directory, self.state_file_name)
            with open(self.state_file_path, 'w') as f:
                f.write(mm.XmlSerializer.serialize(self.state))

            simulation.context.setTime(0)
        else:
            lastpositions = simulation.context.getState(getPositions=True, getVelocities=True, enforcePeriodicBox=False,
                                                        getForces=True, getEnergy=True, getIntegratorParameters=True,
                                                        getParameterDerivatives=True).getPositions()
            last_pdb_file_path = os.path.join(self.output_directory, self.last_pdb_filename)
            self.last_pdb = app.PDBFile.writeFile(modeller.topology, lastpositions, open(last_pdb_file_path, 'w'),
                                                  keepIds=True)
            print("Simulation Just Minimized ...")

            self.state = simulation.context.getState(getPositions=True, getVelocities=True, enforcePeriodicBox=False,
                                                     getForces=True, getEnergy=True,
                                                     getIntegratorParameters=True, getParameterDerivatives=True)

            self.state_file_path = os.path.join(self.output_directory, self.state_file_name)
            with open(self.state_file_path, 'w') as f:
                f.write(mm.XmlSerializer.serialize(self.state))

# Reference_MD_Engine('2j0w.pdb', minimize=True, minimize_steps=50, platform_name='OpenCL', water_padding=10)
