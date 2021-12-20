from functools import partial

import mdtraj
from openmm import app, Platform, VerletIntegrator, VariableVerletIntegrator, LangevinIntegrator, \
    BrownianIntegrator, CustomIntegrator, XmlSerializer, NonbondedForce, MTSIntegrator, CustomNonbondedForce, \
    HarmonicBondForce, HarmonicAngleForce, PeriodicTorsionForce, CustomBondForce, CMMotionRemover, CMAPTorsionForce, \
    CustomTorsionForce

import openmm.unit as u
from sys import stdout
import parmed as pmd
import numpy as np
import time
import os
import pandas as pd
# from ui_main import *
# from src.builder import *

# from .get_positions_from_trajectory_file import get_openmm_pos_from_traj
# from .response_time_creator import getResidueResponseTimes

VELUNIT = u.angstrom / u.picosecond
RESIDUE_NAME = str
CpuThreads = 4


def get_snapshot_positions(pdb_files, pdb_positions=[]):
    for pdb in pdb_files:
        pdb_positions.append(app.PDBFile(pdb).positions)
    return pdb_positions


def collect_all_parameters(system):
    forces = {system.getForce(index).__class__.__name__: system.getForce(
        index) for index in range(system.getNumForces())}

    NonBonded_Parameters = {}
    NonBonded_Exception_Parameters = {}

    Harmonic_Bond_Force_Parameters = {}
    Harmonic_Angle_Force_Parameters = {}
    Periodic_Torsion_Force_Parameters = {}

    Nonbonded_F = forces['NonbondedForce']

    # Harmonic_BF = forces['HarmonicBondForce']
    # Harmonic_AF = forces['HarmonicAngleForce']
    # Periodic_TF = forces['PeriodicTorsionForce']

    for index in range(Nonbonded_F.getNumParticles()):
        charge, sigma, epsilon = Nonbonded_F.getParticleParameters(index)
        NonBonded_Parameters[index] = (charge, sigma, epsilon)

    for index in range(Nonbonded_F.getNumExceptions()):
        particle1, particle2, chargeProd, sigma, epsilon = Nonbonded_F.getExceptionParameters(index)
        NonBonded_Exception_Parameters[index] = (particle1, particle2, chargeProd, sigma, epsilon)

    # for index in range(Harmonic_BF.getNumBonds()):
    #     p1, p2, length, k = Harmonic_BF.getBondParameters(index)
    #     Harmonic_Bond_Force_Parameters[index] = (p1, p2, length, k)
    #
    # for index in range(Harmonic_AF.getNumAngles()):
    #     p1, p2, p3, angle, k = Harmonic_AF.getAngleParameters(index)
    #     Harmonic_Angle_Force_Parameters[index] = (p1, p2, p3, angle, k)
    #
    # for index in range(Periodic_TF.getNumTorsions()):
    #     p1, p2, p3, p4, periodicity, phase, k = Periodic_TF.getTorsionParameters(index)
    #     Periodic_Torsion_Force_Parameters[index] = (p1, p2, p3, p4, periodicity, phase, k)

    # return NonBonded_Parameters, Harmonic_Bond_Force_Parameters, Harmonic_Angle_Force_Parameters, \
    #        Periodic_Torsion_Force_Parameters, Nonbonded_F, Harmonic_BF, Harmonic_AF, Periodic_TF, \
    #        NonBonded_Exception_Parameters
    return NonBonded_Parameters, Harmonic_Bond_Force_Parameters, Harmonic_Angle_Force_Parameters, \
           Periodic_Torsion_Force_Parameters, Nonbonded_F, None, None, None, \
           NonBonded_Exception_Parameters


def get_residue_atoms():
    RESIDUE_NAME = {}
    per_residue_atoms = {}

    residues = [r for r in modeller.topology.residues()]  # build a list of residues

    for index in range(len(residues)):
        indicated_res_atoms = residues[index].atoms()
        RESIDUE_NAME[index] = residues[index].name
        per_residue_atoms[index] = list(indicated_res_atoms)

    return per_residue_atoms, RESIDUE_NAME, len(residues)


def residue_based_decomposition(topol, trj_pos_list, start_res, stop_res, output_directory, ref_energy_name,
                                modif_energy_name, origin_last_pdb, ff, que=None):
    global modeller, by_pass_bonds_index, simulation, by_pass_atoms, nonbonded_group_num, modified_df, reference_df

    print('Loading...')

    protein_ff = ff
    pdb = app.PDBFile(topol)

    modeller = app.Modeller(pdb.topology, pdb.positions)

    residue_length_of_protein = len([r for r in modeller.topology.residues()])
    print("residue length of protein is %s" % residue_length_of_protein)

    if ref_energy_name is None:
        modified_df = pd.DataFrame(
            np.random.randint(1, size=(int(len(trj_pos_list)), residue_length_of_protein))).astype(
            float)

    if ref_energy_name is not None:
        reference_df = pd.DataFrame(
            np.random.randint(1, size=(int(len(trj_pos_list) / 2), residue_length_of_protein))).astype(float)

        modified_df = pd.DataFrame(
            np.random.randint(1, size=(int(len(trj_pos_list) / 2), residue_length_of_protein))).astype(
            float)

    structure = pmd.load_file(origin_last_pdb)
    forcefield = app.ForceField(protein_ff)
    # forcefield = app.ForceField('amber96.xml')

    # state = XmlSerializer.deserialize(open('state.xml').read())
    # box_vectors = Quantity(np.diag([80, 80, 80]), angstrom)
    # modeller.topology.setPeriodicBoxVectors(box_vectors)
    modeller.topology.setPeriodicBoxVectors(structure.box_vectors)
    # print(state.getPeriodicBoxVectors())
    # system = forcefield.createSystem(modeller.topology, nonbondedMethod=app.PME,
    #                                  nonbondedCutoff=12.0 * u.angstrom,
    #                                  ewaldErrorTolerance=1e-5,
    #                                  constraints=None, rigidWater=None)
    system = forcefield.createSystem(modeller.topology, nonbondedMethod=app.NoCutoff, constraints=None,
                                     rigidWater=None)

    # system = forcefield.createSystem(modeller.topology, nonbondedMethod=app.NoCutoff,
    #                                  constraints=None, rigidWater=True)

    # platform = Platform.getPlatformByName('CUDA')
    # properties = {'CudaPrecision': 'mixed', 'DeterministicForces': 'true'}

    platform = Platform.getPlatformByName('OpenCL')
    # properties = {'Threads': '8'}
    properties = {'Precision': 'mixed'}
    # properties = {'OpenCLDeviceIndex': '0'}

    # properties = {  # 'Threads': '1', #number of threads for CPU - all definitions must be strings (I think)
    #     'Precision': 'mixed',  # for CUDA or OpenCL, select the precision (single, mixed, or double)
    #     'DeviceIndex': '0',  # selects which GPUs to use - set this to zero if using CUDA_VISIBLE_DEVICES
    #     'DeterministicForces': 'True'  # Makes sure forces with CUDA and PME are deterministic
    # }

    # platform = Platform.getPlatformByName('CPU')
    # properties = {'CpuThreads': '4'}
    nonbonded_group = None
    for i, f in enumerate(system.getForces()):
        print(i, f.__class__.__name__)
        if isinstance(f, NonbondedForce):
            f.setForceGroup(i)

            # f.setIncludeDirectSpace(False)
            nonbonded_group_num = i  # force.getForceGroup()
            # f.setUseSwitchingFunction(use=True)
            # f.setSwitchingDistance(10.0 * u.angstrom)
            # force.setUseSwitchingFunction(use=True)
            # force.setSwitchingDistance(10.0*u.angstrom)
            f.setUseDispersionCorrection(False)
            print(f.getReactionFieldDielectric())
            f.setReactionFieldDielectric(1.0)

        if isinstance(f, HarmonicBondForce):
            print("Harmonic Bond Force deleted")
            system.removeForce(i)
        if isinstance(f, HarmonicAngleForce):
            print("Harmonic Angle Force deleted")
            system.removeForce(i)

        if isinstance(f, PeriodicTorsionForce):
            print("Periodic Torsion Force deleted")
            system.removeForce(i)

        if isinstance(f, CustomTorsionForce):
            print("Custom Torsion Force deleted")
            system.removeForce(i)

        if isinstance(f, CMAPTorsionForce):
            print("CMAP Torsion Force deleted")
            system.removeForce(i)

    # for force in system.getForces():
    #     force.setForceGroup(0)  # all forces default to group 0
    #     if force.__class__.__name__ == 'NonbondedForce':
    #         # force.setReciprocalSpaceForceGroup(1)
    #         force.setForceGroup(2)
    #         force.setIncludeDirectSpace(False)
    #         nonbonded_group_num = 2 #force.getForceGroup()
    #         force.setUseSwitchingFunction(use=True)
    #         force.setSwitchingDistance(10.0 * u.angstrom)
    #         # force.setUseSwitchingFunction(use=True)
    #         # force.setSwitchingDistance(10.0*u.angstrom)
    #         force.setUseDispersionCorrection(False)

    # for force in system.getForces():
    #     if force.__class__.__name__ == 'NonbondedForce':
    #         print("Force Groups: %s" % force.getForceGroup())
    #         if force.getForceGroup() == 2:
    #             force.setIncludeDirectSpace(False)
    #             nonbonded_group_num = force.getForceGroup()
    #             force.setUseSwitchingFunction(use=True)
    #             force.setSwitchingDistance(10.0 * u.angstrom)
    #             # force.setUseSwitchingFunction(use=True)
    #             # force.setSwitchingDistance(10.0*u.angstrom)
    #             force.setUseDispersionCorrection(False)
    #     # print(force, force.getForceGroup())

    integrator = VerletIntegrator(1.0 * u.femtoseconds)
    integrator.setConstraintTolerance(0.000001)

    simulation = app.Simulation(modeller.topology, system, integrator, platform, properties)

    NonBonded_Parameters, Harmonic_Bond_Force_Parameters, \
    Harmonic_Angle_Force_Parameters, Periodic_Torsion_Force_Parameters, \
    Nonbonded_F, Harmonic_BF, Harmonic_AF, Periodic_TF, NonBonded_Exception_Parameters = collect_all_parameters(system)

    per_res_atoms_list, selected_res_name, res_number = get_residue_atoms()

    simulation.context.setPositions(modeller.positions)
    simulation.context.computeVirtualSites()

    for res_num in range(len(per_res_atoms_list)):  # for res_num in range(len(all_residues)):
        by_pass_atoms_index = []
        for i in range(len(per_res_atoms_list[res_num])):
            by_pass_atoms_index.append(per_res_atoms_list[res_num][i].index)

        # --------------------------------------------- NOBONDED FORCE ----------------------------------------------- #
        if simulation.context.getPlatform().getName() in ['OpenCL', 'CUDA']:
            for index in range(Nonbonded_F.getNumParticles()):
                charge, sigma, epsilon = NonBonded_Parameters[index]
                exclude = (index not in by_pass_atoms_index)
                Nonbonded_F.setParticleParameters(index, 0 * charge if exclude else charge, sigma,
                                                  0 * epsilon if exclude else epsilon)

            for index in range(Nonbonded_F.getNumExceptions()):
                exclude = (index not in by_pass_atoms_index)
                particle1, particle2, chargeProd, sigma, epsilon = NonBonded_Exception_Parameters[index]
                Nonbonded_F.setExceptionParameters(index, particle1, particle2,
                                                   0 * chargeProd if exclude else chargeProd,
                                                   sigma, 0 * epsilon if exclude else epsilon)

            Nonbonded_F.updateParametersInContext(simulation.context)

        if simulation.context.getPlatform().getName() == 'CPU':
            for index in range(Nonbonded_F.getNumParticles()):
                charge, sigma, epsilon = NonBonded_Parameters[index]
                exclude = (index not in by_pass_atoms_index)
                Nonbonded_F.setParticleParameters(index, 0 * charge if exclude else charge, sigma,
                                                  0 * epsilon if exclude else epsilon)

            for index in range(Nonbonded_F.getNumExceptions()):
                exclude = (index not in by_pass_atoms_index)
                particle1, particle2, chargeProd, sigma, epsilon = NonBonded_Exception_Parameters[index]
                Nonbonded_F.setExceptionParameters(index, particle1, particle2,
                                                   10e-15 * chargeProd if exclude else chargeProd,
                                                   sigma, 0 * epsilon if exclude else epsilon)

            Nonbonded_F.updateParametersInContext(simulation.context)

        # ------------------------------------------- HARMONIC BOND FORCE -------------------------------------------- #
        """
        for index in range(Harmonic_BF.getNumBonds()):
            p1, p2, length, k = Harmonic_Bond_Force_Parameters[index]
            exclude = (p1 not in by_pass_atoms_index or p2 not in by_pass_atoms_index)
            Harmonic_BF.setBondParameters(index, p1, p2, length, 0 if exclude else k)

        Harmonic_BF.updateParametersInContext(simulation.context)

        # ------------------------------------------ HARMONIC ANGLES FORCE ------------------------------------------- #

        for index in range(Harmonic_AF.getNumAngles()):
            p1, p2, p3, angle, k = Harmonic_Angle_Force_Parameters[index]
            exclude = (p1 not in by_pass_atoms_index or p2 not in by_pass_atoms_index or p3 not in by_pass_atoms_index)
            Harmonic_AF.setAngleParameters(index, p1, p2, p3, angle, 0 if exclude else k)

        Harmonic_AF.updateParametersInContext(simulation.context)

        # ------------------------------------------ PERIODIC TORSION FORCE ------------------------------------------ #

        for index in range(Periodic_TF.getNumTorsions()):
            p1, p2, p3, p4, periodicity, phase, k = Periodic_Torsion_Force_Parameters[index]
            exclude = (
                    p1 not in by_pass_atoms_index or p2 not in by_pass_atoms_index or p3 not in by_pass_atoms_index or p4 not in by_pass_atoms_index)
            Periodic_TF.setTorsionParameters(index, p1, p2, p3, p4, periodicity, phase, 0 if exclude else k)

        Periodic_TF.updateParametersInContext(simulation.context)
        """

        for i in range(len(trj_pos_list)):
            simulation.context.setPositions(trj_pos_list[i])

            st = simulation.context.getState(getEnergy=True,
                                             groups=1 << nonbonded_group_num).getPotentialEnergy().value_in_unit(
                u.kilocalories_per_mole)

            if ref_energy_name is None:
                modified_df.loc[i][res_num] = float(st)

            if ref_energy_name is not None:
                if i < int(len(trj_pos_list) / 2):
                    reference_df.loc[i][res_num] = float(st)

                else:
                    modified_df.loc[i - int(len(trj_pos_list) / 2)][res_num] = float(st)
        if que is not None:
            que.put([res_num, res_number, st])  # res_num = current residue  -  res_number = all residues number

        # print("#################### %s #####################" % res_num)

        del by_pass_atoms_index

    if ref_energy_name is None:
        modified_df.to_csv(os.path.join(output_directory, modif_energy_name), index=False)

    if ref_energy_name is not None:
        reference_df.to_csv(os.path.join(output_directory, ref_energy_name), index=False)
        modified_df.to_csv(os.path.join(output_directory, modif_energy_name), index=False)

    if que is not None:
        que.put("Progress Finished Succesfully :)")  # res_num = current residue  -  res_number = all residues number
# from get_positions_from_trajectory_file import get_openmm_pos_from_traj_with_mdtraj
#
# position_list, unwrap_pdb = get_openmm_pos_from_traj_with_mdtraj(
#     'C:/Users/HIbrahim/Desktop/perturbation/per_residue_energy/example/last.pdb',
#     'C:/Users/HIbrahim/Desktop/perturbation/per_residue_energy/example/energy_perturbation_trajectory3.dcd',
#     'C:/Users/HIbrahim/Desktop/perturbation/per_residue_energy/example/without_energy_perturbation_trajectory.dcd')
#
# # --> RESIDUE BASED DECOMPOSITION
# residue_based_decomposition(topol=unwrap_pdb, trj_pos_list=position_list, start_res=0, stop_res=250, ff='charmm36',
#                             output_directory=os.getcwd(),
#                             ref_energy_name='reference_energy_file_1.csv',
#                             modif_energy_name='modified_energy_file_1.csv',
#                             origin_last_pdb='C:/Users/HIbrahim/Desktop/perturbation/per_residue_energy/example/last.pdb')
