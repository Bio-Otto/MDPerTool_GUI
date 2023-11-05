from pdbfixer import PDBFixer
from openmm.app.pdbfile import PDBFile
from openmm import app
import openmm.app.data
import openmm.app.data.charmm36
from sys import stdout
import parmed as pmd
import numpy as np
import time
import os


def fix_pdb(pdb_id, fixed_pdb_out_path, pH=7.4, logger_object=None):

    if len(pdb_id) != 4 and logger_object is not None:
        logger_object.info("Creating PDBFixer...".format())
        fixer = PDBFixer(pdb_id)
        logger_object.info("Finding missing residues...".format())
        fixer.findMissingResidues()

        chains = list(fixer.topology.chains())
        keys = fixer.missingResidues.keys()
        for key in list(keys):
            chain = chains[key[0]]
            if key[1] == 0 or key[1] == len(list(chain.residues())):
                del fixer.missingResidues[key]

        logger_object.info("Finding nonstandard residues...".format())
        fixer.findNonstandardResidues()

        logger_object.info("Replacing nonstandard residues...".format())
        fixer.replaceNonstandardResidues()

        logger_object.info("Removing heterogens...".format())
        fixer.removeHeterogens(keepWater=False)

        logger_object.info("Finding missing atoms...".format())
        logger_object.info("Founded missing  Residues: %s" % fixer.missingResidues)
        fixer.findMissingAtoms()

        logger_object.info("Adding missing atoms...".format())
        fixer.addMissingAtoms()
        logger_object.info("Adding missing hydrogens...".format())
        #fixer.addMissingHydrogens(pH)

        logger_object.info("Writing PDB file...".format())
        fixed_pdb_out_path = os.path.join(fixed_pdb_out_path, "%s_fixed_pH_%s.pdb" % (os.path.basename(pdb_id).split('.')[0], pH))
        PDBFile.writeFile(fixer.topology, fixer.positions, open(fixed_pdb_out_path, "w"), keepIds=True)

        return fixed_pdb_out_path

    if len(pdb_id) != 4 and logger_object is None:
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
        print(fixer.missingResidues)
        print("Finding missing atoms...")
        fixer.findMissingAtoms()

        print("Adding missing atoms...")
        fixer.addMissingAtoms()
        print("Adding missing hydrogens...")
        #fixer.addMissingHydrogens(pH)

        print("Writing PDB file...")
        fixed_pdb_out_path = os.path.join(fixed_pdb_out_path, "%s_fixed_pH_%s.pdb" % (os.path.basename(pdb_id).split('.')[0], pH))
        PDBFile.writeFile(fixer.topology, fixer.positions, open(fixed_pdb_out_path, "w"), keepIds=True)

        return fixed_pdb_out_path
