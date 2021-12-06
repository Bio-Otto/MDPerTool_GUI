import pandas as pd
from simtk.openmm import app
from simtk.openmm import *
from simtk.unit import *

from sys import stdout
import parmed as pmd
import numpy as np
from numpy.random import randint


# reference_pdb_traj = 'traj_x4_reference.pdb'
# modified_pdb_traj = 'traj_x4_modified_velocity.pdb'


def read_model_from_multi_pdb(reference_pdb_traj, modified_pdb_traj):
    reference_snapshots_list = []
    modified_snapshots_list = []
    with open(reference_pdb_traj) as reference_multi_pdb:
        reference_multi_pdb_text = reference_multi_pdb.read()

    with open(modified_pdb_traj) as modified_multi_pdb:
        modified_multi_pdb_text = modified_multi_pdb.read()

    reference_model_number = 1
    reference_new_file_text = ""
    for line in filter(None, reference_multi_pdb_text.splitlines()):
        line = line.strip()  # for better control of ends of lines
        if line == "ENDMDL":
            # save file with file number in name
            output_file = open('reference_' + str(reference_model_number) + ".pdb", "w")
            reference_snapshots_list.append('reference_' + str(reference_model_number) + ".pdb")
            output_file.write(reference_new_file_text.rstrip('\r\n'))  # rstrip to remove trailing newline
            output_file.close()
            # reset everything for next model
            reference_model_number += 1
            reference_new_file_text = ""
        elif not (line.startswith("MODEL") or line.startswith("HETATM") or line.startswith("TER")):
            reference_new_file_text += line + '\n'

    modified_model_number = 1
    modified_new_file_text = ""
    for line in filter(None, modified_multi_pdb_text.splitlines()):
        line = line.strip()  # for better control of ends of lines
        if line == "ENDMDL":
            # save file with file number in name
            output_file = open('modified_' + str(modified_model_number) + ".pdb", "w")
            modified_snapshots_list.append('modified_' + str(modified_model_number) + ".pdb")
            output_file.write(modified_new_file_text.rstrip('\r\n'))  # rstrip to remove trailing newline
            output_file.close()
            # reset everything for next model
            modified_model_number += 1
            modified_new_file_text = ""
        elif not (line.startswith("MODEL") or line.startswith("HETATM") or line.startswith("TER")):
            modified_new_file_text += line + '\n'

    return reference_snapshots_list, modified_snapshots_list

# read_model_from_multi_pdb(reference_pdb_traj, modified_pdb_traj)
