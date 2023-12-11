# import MDAnalysis as md
# from openmm import Vec3
# from openmm.unit import nanometer
import mdtraj as md
import os
from sys import stdout
import numpy as np
from pathlib import Path
import sys
import os
import MDAnalysis as mda

"""
def get_openmm_pos_from_traj(top, ref_traj, modif_traj, selected_atoms='protein', write_dcd=True,
                             write_pdb_for_selection=True, start=None, stop=None):
    global first_snapshot_name
    positions = []
    All_Positions = {}

    save_directory = os.path.dirname(ref_traj)

    # set the universe object
    u = md.Universe(top, ref_traj, modif_traj)

    # set atom selection for all atoms of interest
    wanted_atoms = u.select_atoms(selected_atoms)

    # Writer object that spans the whole trajectory and
    # writes an abbreviated trajectory only including atoms described by atom selection
    for cnt in u.trajectory:
        positions.append(wanted_atoms.positions)

    for k, i in enumerate(positions):
        Pos = []
        for xyzi in i:
            Pos.append(Vec3(xyzi[0], xyzi[1], xyzi[2]))

        All_Positions[k] = (Pos * nanometer)
        del Pos

    if write_dcd:
        if start and stop is not None:
            with md.Writer(os.path.join(save_directory, '%s_%s_%s.dcd' % (ref_traj, start, stop)),
                                   len(wanted_atoms.atoms)) as V:
                for ts in u.trajectory[start:stop]:
                    V.write(wanted_atoms)

        if start is None and stop is not None:
            start = 0
            with md.Writer(os.path.join(save_directory, '%s_%s_%s.dcd' % (ref_traj, start, stop)),
                                   len(wanted_atoms.atoms)) as V:
                for ts in u.trajectory[start:stop]:
                    V.write(wanted_atoms)

        if start is not None and stop is None:
            with md.Writer(os.path.join(save_directory, '%s_%s_all.dcd' % (ref_traj, start)),
                                   len(wanted_atoms.atoms)) as V:
                for ts in u.trajectory[start:]:
                    V.write(wanted_atoms)

    if write_pdb_for_selection:
        first_snapshot_name = os.path.join(save_directory, '%s_first_snapshot.pdb' % selected_atoms)
        wanted_atoms.write(first_snapshot_name)

    return All_Positions, first_snapshot_name
"""


def create_restrained_atoms_topology(pdb_file):
    print("Creating topology with only the restrained atoms.")

    file = md.load(pdb_file).remove_solvent()
    topology = file.topology

    box_dimension = file.openmm_boxes(0)

    all_atoms_indices = []

    for i in range(topology.n_residues):
        all_atoms_indices.append(topology.select("resid %s" % i))

    return all_atoms_indices, topology.to_openmm(), box_dimension


def get_openmm_pos_from_traj_with_mdtraj(top, ref_traj, modif_traj, selected_atoms='protein',
                                         write_pdb_for_selection=True, start=None, stop=None, stride=1,
                                         logger_object=None):
    global first_snapshot_name
    positions = []
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

    # set atom selection for all atoms of interest
    wanted_atoms_traj = trajectory_collector.atom_slice(trajectory_collector.topology.select(selected_atoms))

    if logger_object is not None:
        logger_object.info("Started to convert Traj positions for OpenMM.".format())

    for i in range(wanted_atoms_traj.n_frames):
        positions.append(wanted_atoms_traj.openmm_positions(i))

    if write_pdb_for_selection:
        first_snapshot_name = os.path.join(save_directory, '%s_first_snapshot.pdb' % selected_atoms)
        wanted_atoms_traj[0].save_pdb(first_snapshot_name)

    del traj1, traj2, trajectory_collector
    return positions, first_snapshot_name


def get_openmm_pos_from_traj_with_mdanalysis(top, ref_traj, modif_traj, selected_atoms='protein',
                                             write_pdb_for_selection=True, start=None, stop=None, stride=1,
                                             logger_object=None):
    global first_snapshot_name
    positions = []
    save_directory = os.path.dirname(modif_traj)

    trajectory_collector = None
    traj1 = None
    traj2 = None
    if ref_traj is not None:
        traj1 = mda.Universe(top, ref_traj)
        traj2 = mda.Universe(top, modif_traj)
        trajectory_collector = mda.Merge(traj1, traj2)


    if ref_traj is None:
        trajectory_collector = mda.Universe(top, modif_traj)

    # set atom selection for all atoms of interest
    wanted_atoms_traj = trajectory_collector.select_atoms(selected_atoms)

    if logger_object is not None:
        logger_object.info("Started to convert Traj positions for OpenMM.")

    for ts in trajectory_collector.trajectory[start:stop:stride]:
        positions.append(wanted_atoms_traj.positions.copy())

    if write_pdb_for_selection:
        first_snapshot_name = os.path.join(save_directory, f'{selected_atoms}_first_snapshot.pdb')
        wanted_atoms_traj.atoms.write(first_snapshot_name, format='PDB')

    # release resources by deleting Universe and Merge objects
    del traj1, traj2, trajectory_collector

    return positions, first_snapshot_name


"""
traj = 'C:/Users/law5_/Desktop/MDPerTool_GUI/mdpertool/output/without_energy_perturbation_trajectory.dcd'
modif_traj = 'C:/Users/law5_/Desktop/MDPerTool_GUI/mdpertool/output/energy_perturbation_trajectory4.dcd'
pdb_file = 'C:/Users/law5_/Desktop/MDPerTool_GUI/mdpertool/output/last.pdb'
#pos_mdtraj, snap_mdtraj = get_openmm_pos_from_traj_with_mdtraj(pdb_file, traj, modif_traj)
#pos_mdanal, snap_mdanal = get_openmm_pos_from_traj_with_mdanalysis(pdb_file, traj, modif_traj)
"""
# all_atoms_indices, new_topology, box_vectors = create_restrained_atoms_topology(pdb_file)
#
# # # print(new_topology)
# All_positions, pdb = get_openmm_pos_from_traj(pdb_file, traj, write_dcd=False)



#
# file = app.PDBFile(pdb)
# top = file.topology
# pos = file.positions
#
#
#
#
#
# print('Forcefield parameters loading to the simulation system ...')
# forcefield = app.ForceField('charmm36.xml')
# #
# #
# #
# #
# print('Constructing an OpenMM System')
# system = forcefield.createSystem(top, nonbondedMethod=app.PME,
#                                  nonbondedCutoff=12.0 * unit.angstroms, constraints=None,
#                                  rigidWater=False,
#                                  ewaldErrorTolerance=5e-5)
#


# print(table)
# topology = topology.subset(self._reporter.analysis_particle_indices)
# # Use the receptor as an anchor molecule and image the ligand.
# anchor_molecules = [{a for a in topology.atoms if a.index in set(topography.receptor_atoms)}]
# imaged_molecules = [{a for a in topology.atoms if a.index in set(topography.ligand_atoms)}]
# # Initialize trajectory object needed for imaging molecules.
# trajectory = mdtraj.Trajectory(xyz=np.zeros((topology.n_atoms, 3)), topology=topology)
