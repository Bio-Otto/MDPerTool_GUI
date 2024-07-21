from pdbfixer import PDBFixer
from openmm.app.pdbfile import PDBFile
from openmm import app
import openmm.app.data
from openmm.unit import *
from sys import stdout
import parmed as pmd
import numpy as np
import time
import os


def apply_mutation(pdb_path, mut_region=None, chain_id=None):
    """
    Make a mutant protein easy.

        Parameters
        ----------

        pdb_path: Give your pdb whole path to this parameter

        mut_region : list of strings
            Each string must include the resName (original), index,
            and resName (target).  For example, ALA-133-GLY will mutate
            alanine 133 to glycine.

        chain_id : str
            Chain ID to apply mutation.

        Example
        ----------

        mutate('C:/Users/HIbrahim/Desktop/MolDynAnalyze/test/last.pdb', mut_region=['ASP-306-ARG'], chain_id='A')

    """

    try:
        pdb_name = os.path.basename(pdb_path).split('.')[0]
        pdb_directory = os.path.dirname(pdb_path)

        mut_file_name = pdb_name + '_chain' + chain_id + '_' + str(mut_region[0]) + '.pdb'
        mut_file_path = os.path.join(pdb_directory, mut_file_name)

        fixer = PDBFixer(pdb_path)
        fixer.applyMutations(mut_region, chain_id)
        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()

        with open(mut_file_path, 'w') as w_file:
            app.PDBFile.writeFile(fixer.topology, fixer.positions, w_file, keepIds=True)

        return mut_file_path

    except ValueError as error:
        print(error)
        print('Please Check Your Input Parameters !!')


def fix_pdb(pdb_id, pH=7.4, output=None, log_obj=None, mutate=None, mut_region=None, chain_id=None):
    if output is not None:
        path = output
    else:
        path = os.getcwd()

    if len(pdb_id) != 4 and log_obj is not None:

        if mutate:
            mut_file_path = apply_mutation(pdb_path=pdb_id, mut_region=mut_region, chain_id=chain_id)

            log_obj.info("Creating PDBFixer...".format())
            fixer = PDBFixer(mut_file_path)
            log_obj.info("Finding missing residues...".format())
            fixer.findMissingResidues()
            log_obj.info("Finding nonstandard residues...".format())
            fixer.findNonstandardResidues()
            log_obj.info("Replacing nonstandard residues...".format())
            fixer.replaceNonstandardResidues()
            log_obj.info("Removing heterogens...".format())
            fixer.removeHeterogens(keepWater=False)
            log_obj.info("Finding missing atoms...".format())
            fixer.findMissingAtoms()
            log_obj.info("Adding missing atoms...".format())
            fixer.addMissingAtoms()
            log_obj.info("Adding missing hydrogens...".format())
            fixer.addMissingHydrogens(pH)
            log_obj.info("Writing PDB file...".format())

            PDBFile.writeFile(
                fixer.topology,
                fixer.positions,
                open(os.path.join(path, "%s_fixed_pH_%s.pdb" % (mut_file_path.split('.')[0], pH)), "w"), keepIds=True)

            return "%s_fixed_pH_%s.pdb" % (mut_file_path.split('.')[0], pH)

        if not mutate:
            log_obj.info("Creating PDBFixer...".format())
            fixer = PDBFixer(pdb_id)
            log_obj.info("Finding missing residues...".format())
            fixer.findMissingResidues()
            log_obj.info("Finding nonstandard residues...".format())
            fixer.findNonstandardResidues()
            log_obj.info("Replacing nonstandard residues...".format())
            fixer.replaceNonstandardResidues()
            log_obj.info("Removing heterogens...".format())
            fixer.removeHeterogens(keepWater=False)
            log_obj.info("Finding missing atoms...".format())
            fixer.findMissingAtoms()
            log_obj.info("Adding missing atoms...".format())
            fixer.addMissingAtoms()
            log_obj.info("Adding missing hydrogens...".format())
            fixer.addMissingHydrogens(pH)
            log_obj.info("Writing PDB file...".format())

            PDBFile.writeFile(
                fixer.topology,
                fixer.positions,
                open(os.path.join(path, "%s_fixed_pH_%s.pdb" % (pdb_id.split('.')[0], pH)), "w"), keepIds=True)

            return "%s_fixed_pH_%s.pdb" % (pdb_id.split('.')[0], pH)

    if len(pdb_id) != 4 and log_obj is None:

        if mutate:
            mut_file_path = apply_mutation(pdb_path=pdb_id, mut_region=mut_region, chain_id=chain_id)

            log_obj.info("Creating PDBFixer...".format())
            fixer = PDBFixer(mut_file_path)
            log_obj.info("Finding missing residues...".format())
            fixer.findMissingResidues()
            log_obj.info("Finding nonstandard residues...".format())
            fixer.findNonstandardResidues()
            log_obj.info("Replacing nonstandard residues...".format())
            fixer.replaceNonstandardResidues()
            log_obj.info("Removing heterogens...".format())
            fixer.removeHeterogens(keepWater=False)
            log_obj.info("Finding missing atoms...".format())
            fixer.findMissingAtoms()
            log_obj.info("Adding missing atoms...".format())
            fixer.addMissingAtoms()
            log_obj.info("Adding missing hydrogens...".format())
            fixer.addMissingHydrogens(pH)
            log_obj.info("Writing PDB file...".format())

            PDBFile.writeFile(
                fixer.topology,
                fixer.positions,
                open(os.path.join(path, "%s_fixed_pH_%s.pdb" % (mut_file_path.split('.')[0], pH)), "w"), keepIds=True)

            return "%s_fixed_pH_%s.pdb" % (mut_file_path.split('.')[0], pH)

        if not mutate:

            log_obj.info("Creating PDBFixer...".format())
            fixer = PDBFixer(pdb_id)
            log_obj.info("Finding missing residues...".format())
            fixer.findMissingResidues()
            log_obj.info("Finding nonstandard residues...".format())
            fixer.findNonstandardResidues()
            log_obj.info("Replacing nonstandard residues...".format())
            fixer.replaceNonstandardResidues()
            log_obj.info("Removing heterogens...".format())
            fixer.removeHeterogens(keepWater=False)
            log_obj.info("Finding missing atoms...".format())
            fixer.findMissingAtoms()
            log_obj.info("Adding missing atoms...".format())
            fixer.addMissingAtoms()
            log_obj.info("Adding missing hydrogens...".format())
            fixer.addMissingHydrogens(pH)
            log_obj.info("Writing PDB file...".format())

            PDBFile.writeFile(
                fixer.topology,
                fixer.positions,
                open(os.path.join(path, "%s_fixed_pH_%s.pdb" % (pdb_id.split('.')[0], pH)), "w"), keepIds=True)

            return "%s_fixed_pH_%s.pdb" % (pdb_id.split('.')[0], pH)