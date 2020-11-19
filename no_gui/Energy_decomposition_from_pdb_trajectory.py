from simtk.openmm import app, Platform, VerletIntegrator, VariableVerletIntegrator, LangevinIntegrator, \
    BrownianIntegrator, CustomIntegrator, XmlSerializer, NonbondedForce, MTSIntegrator
from simtk.unit import *
from sys import stdout
import parmed as pmd
import numpy as np
import time
import pandas as pd
from .response_time_creator import getResidueResponseTimes

VELUNIT = angstrom / picosecond
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
    Harmonic_BF = forces['HarmonicBondForce']
    Harmonic_AF = forces['HarmonicAngleForce']
    Periodic_TF = forces['PeriodicTorsionForce']

    print("Harmoniv BOND FORCES: %s" % Harmonic_BF.getNumBonds())

    for index in range(Nonbonded_F.getNumParticles()):
        charge, sigma, epsilon = Nonbonded_F.getParticleParameters(index)
        NonBonded_Parameters[index] = (charge, sigma, epsilon)

    for index in range(Nonbonded_F.getNumExceptions()):
        particle1, particle2, chargeProd, sigma, epsilon = Nonbonded_F.getExceptionParameters(index)
        NonBonded_Exception_Parameters[index] = (particle1, particle2, chargeProd, sigma, epsilon)

    for index in range(Harmonic_BF.getNumBonds()):
        p1, p2, length, k = Harmonic_BF.getBondParameters(index)
        Harmonic_Bond_Force_Parameters[index] = (p1, p2, length, k)

    for index in range(Harmonic_AF.getNumAngles()):
        p1, p2, p3, angle, k = Harmonic_AF.getAngleParameters(index)
        Harmonic_Angle_Force_Parameters[index] = (p1, p2, p3, angle, k)

    for index in range(Periodic_TF.getNumTorsions()):
        p1, p2, p3, p4, periodicity, phase, k = Periodic_TF.getTorsionParameters(index)
        Periodic_Torsion_Force_Parameters[index] = (p1, p2, p3, p4, periodicity, phase, k)

    return NonBonded_Parameters, Harmonic_Bond_Force_Parameters, Harmonic_Angle_Force_Parameters, \
           Periodic_Torsion_Force_Parameters, Nonbonded_F, Harmonic_BF, Harmonic_AF, Periodic_TF, \
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


def main(pdb_filename, pos_list, start_res, stop_res):
    global modeller, by_pass_bonds_index, simulation, by_pass_atoms, nonbonded_group_num
    print('Loading...')

    pdb = app.PDBFile(pdb_filename)

    modeller = app.Modeller(pdb.topology, pdb.positions)

    print("here")
    residue_length_of_protein = len([r for r in modeller.topology.residues()])
    print("residue length of protein is %s" % residue_length_of_protein)

    reference_df = pd.DataFrame(
        np.random.randint(1, size=(int(len(pos_list) / 2), residue_length_of_protein))).astype(float)

    modified_df = pd.DataFrame(
        np.random.randint(1, size=(int(len(pos_list) / 2), residue_length_of_protein))).astype(
        float)

    structure = pmd.load_file(pdb_filename)
    forcefield = app.ForceField('charmm36.xml')
    # forcefield = app.ForceField('amber96.xml')

    # state = XmlSerializer.deserialize(open('state.xml').read())
    # box_vectors = Quantity(np.diag([80, 80, 80]), angstrom)
    # modeller.topology.setPeriodicBoxVectors(box_vectors)
    modeller.topology.setPeriodicBoxVectors(structure.box_vectors)
    # print(state.getPeriodicBoxVectors())
    # system = forcefield.createSystem(modeller.topology, nonbondedMethod=app.PME, ewaldErrorTolerance=5e-5,
    #                                  constraints=None, rigidWater=False)
    system = forcefield.createSystem(modeller.topology, nonbondedMethod=app.PME, constraints=None,
                                     rigidWater=False)
    print("here")

    platform = Platform.getPlatformByName('OpenCL')
    properties = {'OpenCLPrecision': 'single'}
    # properties = {'OpenCLDeviceIndex': '0'}

    # platform = Platform.getPlatformByName('CPU')
    # properties = {'CpuThreads': '4'}


    for i, f in enumerate(system.getForces()):
        if isinstance(f, NonbondedForce):
            f.setReciprocalSpaceForceGroup(i)
        f.setForceGroup(i)


    for force in system.getForces():
        if force.__class__.__name__ == 'NonbondedForce':
            nonbonded_group_num = force.getForceGroup()
        print(force, force.getForceGroup())

    integrator = VerletIntegrator(0.001)

    simulation = app.Simulation(modeller.topology, system, integrator, platform, properties)

    NonBonded_Parameters, Harmonic_Bond_Force_Parameters, \
    Harmonic_Angle_Force_Parameters, Periodic_Torsion_Force_Parameters, \
    Nonbonded_F, Harmonic_BF, Harmonic_AF, Periodic_TF, NonBonded_Exception_Parameters = collect_all_parameters(system)

    # log = open('log_%s.txt' % model_number, 'a')
    per_res_atoms_list, selected_res_name, res_number = get_residue_atoms()

    # per_res_atoms_list, selected_res_name = get_residue_atoms()

    simulation.context.setPositions(modeller.positions)
    simulation.context.computeVirtualSites()

    for res_num in range(res_number):  # for res_num in range(len(all_residues)):

        by_pass_atoms_index = []
        for i in range(len(per_res_atoms_list[res_num])):
            by_pass_atoms_index.append(per_res_atoms_list[res_num][i].index)

        ######################################################################### NOBONDED FORCE
        if simulation.context.getPlatform().getName() == 'OpenCL':
            print("OpenCL Kullanımda")
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
            print("CPU Kullanımda")
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

        ######################################################################### HARMONIC BOND FORCE

        for index in range(Harmonic_BF.getNumBonds()):
            p1, p2, length, k = Harmonic_Bond_Force_Parameters[index]
            exclude = (p1 not in by_pass_atoms_index or p2 not in by_pass_atoms_index)
            Harmonic_BF.setBondParameters(index, p1, p2, length, 0 if exclude else k)

        Harmonic_BF.updateParametersInContext(simulation.context)

        ######################################################################### HARMONIC ANGLES FORCE

        for index in range(Harmonic_AF.getNumAngles()):
            p1, p2, p3, angle, k = Harmonic_Angle_Force_Parameters[index]
            exclude = (p1 not in by_pass_atoms_index or p2 not in by_pass_atoms_index or p3 not in by_pass_atoms_index)
            Harmonic_AF.setAngleParameters(index, p1, p2, p3, angle, 0 if exclude else k)
        Harmonic_AF.updateParametersInContext(simulation.context)
        ######################################################################### HARMONIC ANGLES FORCE

        for index in range(Periodic_TF.getNumTorsions()):
            p1, p2, p3, p4, periodicity, phase, k = Periodic_Torsion_Force_Parameters[index]
            exclude = (
                    p1 not in by_pass_atoms_index or p2 not in by_pass_atoms_index or p3 not in by_pass_atoms_index or p4 not in by_pass_atoms_index)
            Periodic_TF.setTorsionParameters(index, p1, p2, p3, p4, periodicity, phase, 0 if exclude else k)
        Periodic_TF.updateParametersInContext(simulation.context)


        for i in range(len(pos_list)):
            simulation.context.setPositions(pos_list[i])
            # st_time_input = time.time()
            print("############################################")

            st = simulation.context.getState(getEnergy=True, groups=1 << nonbonded_group_num).getPotentialEnergy().value_in_unit(
                kilocalories_per_mole)
            print(st)
            print("RES_NUM: %s" % res_num)


            if i < int(len(pos_list) / 2):
                reference_df.loc[i][res_num] = float(st)

            else:
                modified_df.loc[i - int(len(pos_list) / 2)][res_num] = float(st)

        print("#################### %s #####################" % res_num)

        del by_pass_atoms_index

    reference_df.to_csv('reference_energy_file.csv', index=False)
    modified_df.to_csv('modified_energy_file.csv', index=False)



