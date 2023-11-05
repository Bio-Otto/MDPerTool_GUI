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


def fix_pdb(pdb_id, output=None, log_obj=None):
    if output is not None:
        path = output
    else:
        path = os.getcwd()

    if len(pdb_id) != 4 and log_obj is not None:
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
        fixer.addMissingHydrogens(7)
        log_obj.info("Writing PDB file...".format())

        PDBFile.writeFile(
            fixer.topology,
            fixer.positions,
            open(os.path.join(path, "%s_fixed_pH_%s.pdb" % (pdb_id.split('.')[0], 7)), "w"), keepIds=True)

        return "%s_fixed_pH_%s.pdb" % (pdb_id.split('.')[0], 7)

    if len(pdb_id) != 4 and log_obj is None:
        print("Creating PDBFixer...")
        fixer = PDBFixer(pdb_id)
        print("Finding missing residues...")
        fixer.findMissingResidues()
        print("Finding nonstandard residues...")
        fixer.findNonstandardResidues()
        print("Replacing nonstandard residues...")
        fixer.replaceNonstandardResidues()
        print("Removing heterogens...")
        fixer.removeHeterogens(keepWater=False)
        print("Finding missing atoms...")
        fixer.findMissingAtoms()
        print("Adding missing atoms...")
        fixer.addMissingAtoms()
        print("Adding missing hydrogens...")
        fixer.addMissingHydrogens(7)
        print("Writing PDB file...")

        PDBFile.writeFile(
            fixer.topology,
            fixer.positions,
            open(os.path.join(path, "%s_fixed_pH_%s.pdb" % (pdb_id.split('.')[0], 7)), "w"), keepIds=True)

        return "%s_fixed_pH_%s.pdb" % (pdb_id.split('.')[0], 7)
