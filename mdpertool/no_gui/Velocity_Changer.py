import pandas as pd
import matplotlib.pyplot as plt
import mdtraj as md
import openmm as mm
import xml.etree.ElementTree as ET
import os


def convert_res_to_atoms(pdb_path, selected_res, atom_to_extract=None):
    """
    Convert residue selections to atom indices.

    Parameters:
    -----------
    pdb_path : str
        Path to the PDB file
    selected_res : list
        List of residues in format ['RES46'] where RES is residue name
    atom_to_extract : str or None
        Specific atom name to extract (e.g., 'CA' for alpha carbon)

    Returns:
    --------
    list
        List of atom indices corresponding to the selected residues
    """
    # Ensure selected_res is a list
    if isinstance(selected_res, str):
        selected_res = [selected_res]

    # Load the topology
    print("PDB-PATH: ", pdb_path)
    try:
        traj = md.load(pdb_path)
        topology = traj.topology

        selected_res_atoms = []

        # Create a lookup of residue strings for comparison
        residue_lookup = {str(res): res for res in topology.residues}

        # Find matching residues
        for res_str in selected_res:
            if res_str in residue_lookup:
                residue = residue_lookup[res_str]
                # Get all atoms for this residue
                for atom in residue.atoms:
                    # If atom_to_extract is specified, only add matching atoms
                    if atom_to_extract is None or atom.name == atom_to_extract:
                        selected_res_atoms.append(atom.index)

        return selected_res_atoms

    except Exception as e:
        print(f"Error processing PDB file: {str(e)}")
        raise


def change_velocity(xml_file, r_factor, modify_atoms):
    tree = ET.parse(xml_file)  # Path to input file
    root = tree.getroot()
    for count, type_tag in enumerate(root.findall('Velocities/Velocity')):

        if count + 1 in modify_atoms:
            type_tag.set('x', str(float(type_tag.get('x')) * r_factor).strip())
            type_tag.set('y', str(float(type_tag.get('y')) * r_factor).strip())
            type_tag.set('z', str(float(type_tag.get('z')) * r_factor).strip())

    save_directory = os.path.join(os.path.dirname(xml_file), 'out_x%s.xml' % r_factor)
    tree.write(save_directory, encoding="UTF-8")
    return save_directory


def read_xml_file(filename):
    """Read serialized object from XML file."""
    with open(filename, 'rt') as infile:
        return mm.XmlSerializer.deserialize(infile.read())

# modify_atoms = convert_res_to_atoms(path, 'SER345', 'CA')
# print(modify_atoms)
# name_of_changed_state_xml = change_velocity('state.xml', 4, modify_atoms)
# print(name_of_changed_state_xml)


