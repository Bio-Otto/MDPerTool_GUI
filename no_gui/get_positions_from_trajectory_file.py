
import MDAnalysis
from simtk.openmm import Vec3
from simtk.unit import nanometer
import mdtraj as md

from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit
from sys import stdout
from simtk.unit import *
from simtk.openmm import *
from simtk.openmm import Context
import numpy as np

def get_openmm_pos_from_traj(top, ref_traj, modif_traj, selected_atoms='protein', write_dcd=True, write_pdb_for_selection=True,
                             start=None, stop=None):

    global first_snapshot_name
    positions = []
    All_Positions = {}

    # set the universe object
    u = MDAnalysis.Universe(top, ref_traj, modif_traj)

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
            with MDAnalysis.Writer('%s_%s_%s.dcd' % (ref_traj, start, stop), len(wanted_atoms.atoms)) as V:
                for ts in u.trajectory[start:stop]:
                    V.write(wanted_atoms)

        if start is None and stop is not None:
            start = 0
            with MDAnalysis.Writer('%s_%s_%s.dcd' % (ref_traj, start, stop), len(wanted_atoms.atoms)) as V:

                for ts in u.trajectory[start:stop]:
                    V.write(wanted_atoms)
        if start is not None and stop is None:
            with MDAnalysis.Writer('%s_%s_all.dcd' % (ref_traj, start), len(wanted_atoms.atoms)) as V:

                for ts in u.trajectory[start:]:
                    V.write(wanted_atoms)

    if write_pdb_for_selection:
        first_snapshot_name = '%s_first_snapshot.pdb' % selected_atoms
        wanted_atoms.write(first_snapshot_name)

    return All_Positions, first_snapshot_name


def create_restrained_atoms_topology(pdb_file):
    print("Creating topology with only the restrained atoms.")

    file = md.load(pdb_file).remove_solvent()
    topology = file.topology

    box_dimension = file.openmm_boxes(0)


    all_atoms_indices = []

    for i in range(topology.n_residues):

        all_atoms_indices.append(topology.select("resid %s" % i))

    return all_atoms_indices, topology.to_openmm(), box_dimension




# traj = 'trajectory.dcd'
# pdb_file = 'last_structure.pdb'



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
