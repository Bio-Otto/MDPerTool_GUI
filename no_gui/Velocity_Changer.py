import pandas as pd
import matplotlib.pyplot as plt
import mdtraj as md
from simtk import openmm
import xml.etree.ElementTree as ET

# path = 'last.pdb'


def convert_res_to_atoms(pdb_path, selected_res, atom_to_extract=None):
    selected_res_atoms = []
    topology = md.load(pdb_path).topology

    for inx, residue in enumerate(topology.residues):
        if str(residue) in selected_res:
            for i in range(len([atom for atom in topology.atoms])):
                if str(topology.atom(i).residue) == selected_res[selected_res.index(str(residue))]:
                    # print(topology.atom(i).residue)
                    selected_res_atoms.append(i)

    if atom_to_extract is not None:
        for i in selected_res_atoms:
            if topology.atom(i - 1).name == atom_to_extract:
                print(i)

    # for i in selected_res_atoms:
    #     print(topology.atom(i).residue)

    return selected_res_atoms


def change_velocity(xml_file, r_factor, modify_atoms):
    tree = ET.parse(xml_file)  # Path to input file
    root = tree.getroot()
    for count, type_tag in enumerate(root.findall('Velocities/Velocity')):

        if count + 1 in modify_atoms:
            type_tag.set('x', str(float(type_tag.get('x')) * r_factor).strip())
            type_tag.set('y', str(float(type_tag.get('y')) * r_factor).strip())
            type_tag.set('z', str(float(type_tag.get('z')) * r_factor).strip())


    tree.write('out_x%s.xml' % r_factor, encoding="UTF-8")
    return 'out_x%s.xml' % r_factor


def read_xml_file(filename):
    """Read serialized object from XML file."""
    with open(filename, 'rt') as infile:
        return openmm.XmlSerializer.deserialize(infile.read())

# modify_atoms = convert_res_to_atoms(path, 'SER345', 'CA')
# print(modify_atoms)
# name_of_changed_state_xml = change_velocity('state.xml', 4, modify_atoms)
# print(name_of_changed_state_xml)


