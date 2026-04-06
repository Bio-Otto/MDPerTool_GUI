import pandas as pd
import matplotlib.pyplot as plt
import mdtraj as md
import openmm as mm
import xml.etree.ElementTree as ET
import os
import re


def _parse_selected_residue_token(token):
    """Parse residue token like ASN17A, ASN17-A, or ASN17."""
    raw = str(token or '').strip()
    match = re.match(r'^([A-Za-z]{3})(\d+)(?:[-_]?([A-Za-z]))?$', raw)
    if not match:
        return None
    return {
        'name': match.group(1).upper(),
        'resseq': int(match.group(2)),
        'chain': (match.group(3) or '').upper(),
    }


def _safe_chain_id_from_residue(residue):
    """Best-effort chain identifier extraction from MDTraj residue object."""
    chain_obj = getattr(residue, 'chain', None)
    if chain_obj is None:
        return ''

    for attr_name in ('chain_id', 'id', 'name'):
        chain_value = getattr(chain_obj, attr_name, None)
        if chain_value is None:
            continue
        chain_text = str(chain_value).strip()
        if chain_text:
            return chain_text.upper()

    chain_index = getattr(chain_obj, 'index', None)
    if chain_index is None:
        return ''
    try:
        # Common convention for simple multimer chains.
        return chr(ord('A') + int(chain_index))
    except Exception:
        return ''


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
        parsed_tokens = [_parse_selected_residue_token(token) for token in selected_res]

        for parsed_token, raw_token in zip(parsed_tokens, selected_res):
            if parsed_token is None:
                # Keep backward compatibility with legacy residue string style.
                for residue in topology.residues:
                    if str(residue).strip() != str(raw_token).strip():
                        continue
                    for atom in residue.atoms:
                        if atom_to_extract is None or atom.name == atom_to_extract:
                            selected_res_atoms.append(atom.index)
                continue

            target_name = parsed_token['name']
            target_resseq = parsed_token['resseq']
            target_chain = parsed_token['chain']

            for residue in topology.residues:
                residue_name = str(getattr(residue, 'name', '')).upper()
                residue_resseq = int(getattr(residue, 'resSeq', -1))
                residue_chain = _safe_chain_id_from_residue(residue)

                if residue_name != target_name or residue_resseq != target_resseq:
                    continue
                if target_chain and residue_chain != target_chain:
                    continue

                for atom in residue.atoms:
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


