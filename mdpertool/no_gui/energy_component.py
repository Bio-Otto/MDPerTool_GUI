import os
import openmm.unit as unit
from openmm.app import CutoffNonPeriodic, Simulation, Modeller, ForceField
from openmm import NonbondedForce, LangevinMiddleIntegrator, Platform
from sys import stdout
import mdtraj as md
import numpy as np
import parmed as pmd
from pdbfixer import PDBFixer
import csv
import pandas as pd


def fix_pdb_file(pdb_id):
    fixer = PDBFixer(pdb_id)
    fixer.removeHeterogens(keepWater=False)

    return fixer


def get_residue_pos_from_traj_with_mdtraj(top, ref_traj, modif_traj, selected_atoms='protein',
                                         write_pdb_for_selection=True, start=None, stop=None, stride=1,
                                         logger_object=None):
    global first_snapshot_name

    All_positions = []
    save_directory = os.path.dirname(modif_traj)

    trajectory_collector = object
    traj1 = object
    traj2 = object
    if ref_traj is not None:
        traj1 = md.join(md.iterload(ref_traj, chunk=100, stride=stride, atom_indices=None, top=top))
        traj2 = md.join(md.iterload(modif_traj, chunk=100, stride=stride, atom_indices=None, top=top))
        trajectory_collector = traj1 + traj2

    if ref_traj is None:
        trajectory_collector = md.join(md.iterload(modif_traj, chunk=100, stride=stride, atom_indices=None, top=top))

    for j in range(selected_atoms):
        positions = []
        # set atom selection for all atoms of interest
        wanted_atoms_traj = trajectory_collector.atom_slice(trajectory_collector.topology.select("resid %s" % j))

        if logger_object is not None:
            logger_object.info("Started to convert Traj positions for OpenMM.".format())

        for i in range(wanted_atoms_traj.n_frames):
            positions.append(wanted_atoms_traj.openmm_positions(i)),

        All_positions.append(positions)
    if write_pdb_for_selection:
        first_snapshot_name = os.path.join(save_directory, '%s_first_snapshot.pdb' % selected_atoms)
        wanted_atoms_traj[0].save_pdb(first_snapshot_name)

    del traj1, traj2, trajectory_collector

    return All_positions, first_snapshot_name, wanted_atoms_traj.n_frames


class EnergyDecomposition:
    def __init__(self, topology_openmm, system, platform, properties, selection_md="protein",
                 nonbondedCutoff=1.2 * unit.nanometers, swithing_distance=1.2 * unit.nanometers, use_switchin_dist=True):
        """The class allows calculating MD energy for a given selection.
            To init the class it's requires to pass openmm topology and system.
            Arguments:
            topology_openmm: Topology of the full system from the openmm
            system: System object from the OpenMM
            selection_md: Selection for the calculated energy, according to the MDTraj syntax
        """
        self.selection_md = selection_md
        # create MDtraj topology
        topology = md.Topology.from_openmm(topology_openmm)
        # select atoms according to the selection
        self.sub_ind = topology.select(self.selection_md)
        # create a slice of topology
        sub_top = topology.subset(self.sub_ind)
        # save old topology to make a from positions
        self.old_topology = topology
        # convert sliced topology to the openmm format
        self.topology = sub_top.to_openmm()

        # Creating system only for protein
        struct = pmd.openmm.load_topology(topology_openmm, system)
        # select structure
        struct = struct[self.sub_ind]
        # if there're hbond restrains, add parameters
        new_bond_type = pmd.topologyobjects.BondType(k=400, req=1.0)
        constrained_bond_type = struct.bond_types.append(new_bond_type)
        struct.bond_types.claim()

        for bond in struct.bonds:
            if bond.type is None:
                bond.type = new_bond_type

        # crete a new system
        new_system = struct.createSystem(nonbondedMethod=CutoffNonPeriodic, nonbondedCutoff=12.0 * pmd.unit.nanometers)

        self.nonbonded_group_num = None
        for i, f in enumerate(new_system.getForces()):
            if isinstance(f, NonbondedForce):
                f.setForceGroup(i)

                # f.setIncludeDirectSpace(False)
                self.nonbonded_group_num = i  # force.getForceGroup()

                if use_switchin_dist:
                    f.setUseSwitchingFunction(use=True)
                    f.setSwitchingDistance(swithing_distance)
                    # force.setUseSwitchingFunction(use=True)
                    # force.setSwitchingDistance(10.0*u.angstrom)
                    # f.setUseDispersionCorrection(False)
                    # f.setReactionFieldDielectric(1.0)
            else:
                f.setForceGroup(i)

        self.system = new_system
        # creating simulation object
        integrator = LangevinMiddleIntegrator(300 * unit.kelvin, 1 / unit.picosecond, 0.004 * unit.picoseconds)
        self.simulation = Simulation(self.topology, self.system, integrator, platform, properties)

    def calc_energy(self, positions, len_of_traj):
        # trajectory = md.Trajectory(positions, self.old_topology)
        # trajectory = md.load(positions, top=native_top)
        # len_of_traj = trajectory.n_frames
        data = []

        # new_positions = [trajectory.atom_slice(self.sub_ind).xyz[i] for i in range(len_of_traj)]
        for i in range(len_of_traj):
            # new_positions = trajectory.atom_slice(self.sub_ind).xyz[i]
            self.simulation.context.setPositions(positions[i])
            # state = self.simulation.context.getState(getEnergy=True)
            # print(state.getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole))

            data.append(self.simulation.context.getState(getEnergy=True, groups=1 << self.nonbonded_group_num).
                        getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole))
            # data.append(self.simulation.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole))

        return data


def process_energy_data(topology_file, protein_ff, water_ff, ref_trajectory, modif_trajectory,
                        platform_type='CPU', nonbonded_cutoff=12.0 * unit.nanometers, ref_energy_name=None,
                        output_directory=None, modif_energy_name=None, device_id_active=True, num_of_threads=1,
                        que=None, logger_object=None):

    if logger_object is not None:
        logger_object.info("The system parameters are being loaded for the decomposition process.".format())



    ############################################################################3

    platform = Platform.getPlatformByName(platform_type)

    if platform_type == 'OpenCL' and device_id_active == True:
        properties = {'OpenCLPrecision': 'double', 'OpenCLDeviceIndex': '1'}
        precision = 'mixed'

    if platform_type == 'OpenCL' and device_id_active == False:
        properties = {'OpenCLPrecision': 'double'}
        precision = 'mixed'

    if platform_type == 'CUDA' and device_id_active == True:
        properties = {'CudaPrecision': 'double', 'CudaDeviceIndex': '1'}
        precision = 'mixed'

    if platform_type == 'CUDA' and device_id_active == False:
        properties = {'CudaPrecision': 'double'}
        precision = 'mixed'

    if platform_type == 'CPU' and device_id_active == True:
        if logger_object is not None:
            logger_object.info("The CPU platform always uses 'mixed' precision.".format())
            logger_object.info("Simulation process will use %s Thread(s)" % num_of_threads)
        properties = {'CpuThreads': '%s' % num_of_threads}
        precision = 'mixed'

    if platform_type == 'Reference':
        if logger_object is not None:
            logger_object.info("The Reference platform always uses 'double' precision.")
        properties = None
        precision = 'double'


    if logger_object is not None:
        logger_object.info("The %s platform will be used for the decomposition process with %s precision."
                           % (platform_type, precision))

    ###########################################################################

    # Fix PDB file
    pdb_top = fix_pdb_file(topology_file)

    # Create a Modeller
    modeller = Modeller(pdb_top.topology, pdb_top.positions)

    # Load structure and force field
    structure = pmd.load_file(topology_file)
    forcefield = ForceField(protein_ff, water_ff)

    # Set periodic box vectors
    modeller.topology.setPeriodicBoxVectors(structure.box_vectors)

    # Create OpenMM system
    system = forcefield.createSystem(modeller.topology, nonbondedMethod=CutoffNonPeriodic, nonbondedCutoff=nonbonded_cutoff)

    # Get residues and create a list of residue names
    residues = modeller.topology.residues()
    res_name_list = [residue.name + residue.id for residue in residues]




    # Get positions from trajectory using custom function (get_openmm_pos_from_traj_with_mdtraj)
    positions, _, len_of_traj = get_residue_pos_from_traj_with_mdtraj(top=topology_file, ref_traj=ref_trajectory,
                                                                     modif_traj=modif_trajectory,
                                                                     selected_atoms=modeller.topology.getNumResidues(),
                                                                     write_pdb_for_selection=True,
                                                                     start=None,
                                                                     stop=None, stride=1,
                                                                     logger_object=None)

    ########################

    if ref_energy_name is None:
        # Dictionary to store energy data for each residue
        modified_df = {res_name: [] for res_name in res_name_list}

    if ref_energy_name is not None:
        # Dictionary to store energy data for each residue
        reference_df = {res_name: [] for res_name in res_name_list}
        # Dictionary to store energy data for each residue
        modified_df = {res_name: [] for res_name in res_name_list}

    ########################


    # Calculate energy decomposition for each residue
    for j in range(modeller.topology.getNumResidues()):
        res_name = res_name_list[j]
        print("Name of Residue: ", res_name)
        enesel = EnergyDecomposition(modeller.topology, system, platform=platform, properties=properties,
                                     selection_md="resid %s" % j, swithing_distance=1.0 * unit.nanometers,
                                     use_switchin_dist=True)
        res_data = enesel.calc_energy(positions[j], len_of_traj)


        #############################
        if ref_energy_name is None:
            modified_df[res_name] = res_data[:int(len_of_traj/2)]

        if ref_energy_name is not None:
            reference_df[res_name] = res_data[:int(len_of_traj / 2)]
            modified_df[res_name] = res_data[int(len_of_traj/2):]
        #############################

        print(res_data)



########################################################################################################################

        if logger_object is not None:
            """
            logger_object.info("Decomposed Res Num: %s, Total Res Num: %s, Decomposition Progress: %s"
                               % (res_num, res_number, res_num/res_number*100))
                               """
            logger_object.log(35, "Decomposed Res Num: %s, Total Res Num: %s, Decomposition Progress: %s"
                              % (j, len(res_name_list), j / len(res_name_list) * 100))
            print("Decomposed Res Num: %s, Total Res Num: %s, Decomposition Progress: %s"
                              % (j, len(res_name_list), j / len(res_name_list) * 100))
########################################################################################################################

    # Write data to a CSV file
    with open(os.path.join(output_directory, ref_energy_name), 'w', newline='') as ref_csv_file:
        writer = csv.writer(ref_csv_file)

        # Write header
        header = list(reference_df.keys())
        writer.writerow(header)

        # Write data to columns
        max_length = max(len(reference_df[key]) for key in header)
        for i in range(max_length):
            row = [reference_df[key][i] if i < len(reference_df[key]) else '' for key in header]
            writer.writerow(row)

        print(f"Data has been written to {ref_energy_name}.")

    # Write data to a CSV file
    with open(os.path.join(output_directory, modif_energy_name), 'w', newline='') as modified_csv_file:
        writer = csv.writer(modified_csv_file)

        # Write header
        header = list(modified_df.keys())
        writer.writerow(header)

        # Write data to columns
        max_length = max(len(modified_df[key]) for key in header)
        for i in range(max_length):
            row = [modified_df[key][i] if i < len(modified_df[key]) else '' for key in header]
            writer.writerow(row)

        print(f"Data has been written to {modif_energy_name}.")

    if logger_object is not None:
        logger_object.info("Progress Finished Succesfully :)")
#
# process_energy_data(topology_file="last.pdb", protein_ff='amber03.xml', water_ff='tip3p.xml',
#                     modif_trajectory='without_energy_perturbation_trajectory.dcd', device_id_active={{Device_ID_active}},
#                     ref_trajectory='energy_perturbation_trajectory4.dcd', platform_type='CPU',
#                     nonbonded_cutoff=1.2 * unit.nanometers, ref_energy_name = 'reference_energy_file_1.csv',
#                     modif_energy_name='modified_energy_file_1.csv', num_of_threads=1, que=None, logger_object=None)

