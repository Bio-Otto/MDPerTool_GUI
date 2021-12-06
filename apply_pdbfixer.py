from pdbfixer import PDBFixer
from simtk.openmm.app.pdbfile import PDBFile
from simtk.openmm import app
import simtk.openmm.app.data
from simtk.unit import *
from sys import stdout
import parmed as pmd
import numpy as np
import time
import os


def fix_pdb(pdb_id, output=None):
    if output is not None:
        path = output
        print("1. de: %s" % path)
    else:
        path = os.getcwd()
        print("2. de: %s" % path)

    if len(pdb_id) != 4:
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
