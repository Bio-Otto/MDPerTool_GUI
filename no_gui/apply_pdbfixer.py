from pdbfixer import PDBFixer
from simtk.openmm.app.pdbfile import PDBFile
from simtk.openmm import app
import simtk.openmm.app.data
import simtk.openmm.app.data.charmm36
from simtk.unit import *
from sys import stdout
import parmed as pmd
import numpy as np
import time
import os


def fix_pdb(pdb_id, fixed_pdb_out_path, pH=7):

    if len(pdb_id) != 4:
        print("Creating PDBFixer...")
        fixer = PDBFixer(pdb_id)
        print("Finding missing residues...")
        fixer.findMissingResidues()

        chains = list(fixer.topology.chains())
        keys = fixer.missingResidues.keys()
        for key in list(keys):
            chain = chains[key[0]]
            if key[1] == 0 or key[1] == len(list(chain.residues())):
                del fixer.missingResidues[key]

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
        fixer.addMissingHydrogens(pH)

        print("Writing PDB file...")
        fixed_pdb_out_path = os.path.join(fixed_pdb_out_path, "%s_fixed_pH_%s.pdb" % (pdb_id.split('.')[0], pH))
        PDBFile.writeFile(fixer.topology, fixer.positions, open(fixed_pdb_out_path, "w"), keepIds=True)

        return fixed_pdb_out_path
