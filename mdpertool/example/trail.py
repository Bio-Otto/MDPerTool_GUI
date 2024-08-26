import numpy as np
import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt
from pdbfixer import PDBFixer

"""
def calculate_magnitude(force):
    x, y, z = force
    return np.sqrt(x ** 2 + y ** 2 + z ** 2)


def read_forces_from_xml(file_path, step):
    forces = []
    try:
        tree = ET.parse(file_path)
        root = tree.getroot()
        for force_element in root.findall('Velocity'):
            if int(force_element.get('step')) == step:
                for atom in force_element.findall('Atom'):
                    x = float(atom.get('x'))
                    y = float(atom.get('y'))
                    z = float(atom.get('z'))
                    forces.append([x, y, z])
                break
    except ET.ParseError as e:
        print(f"XML parsing error: {e}")
    except Exception as e:
        print(f"Error: {e}")
    return np.array(forces)


def track_significant_atoms(ref_file_path, dis_file_path, steps, threshold=2000):
    significant_atoms = {}  # Dictionary to store step -> atom indices

    for step in steps:
        ref_forces = read_forces_from_xml(ref_file_path, step)
        dis_forces = read_forces_from_xml(dis_file_path, step)

        if ref_forces.size > 0 and dis_forces.size > 0 and ref_forces.shape == dis_forces.shape:
            forces_difference = dis_forces - ref_forces
            for i, force_diff in enumerate(forces_difference):
                mag = calculate_magnitude(force_diff)
                if mag >= threshold:
                    if step not in significant_atoms:
                        significant_atoms[step] = []
                    significant_atoms[step].append((i, mag))

    return significant_atoms


def plot_significant_pathway(significant_atoms):
    for step, atoms in significant_atoms.items():
        for atom_index, mag in atoms:
            plt.scatter(step, atom_index, s=mag / 100, label=f'Atom {atom_index}')
    plt.xlabel('Step')
    plt.ylabel('Atom Index')
    plt.title('Significant Atom Pathway')
    plt.legend(loc='best')
    plt.show()


# Example usage:
ref_file_path = 'C:\\Users\\law5_\\Desktop\\MDPerTool_GUI\\mdpertool\\output\\ref_protein_velocities.xml'
dis_file_path = 'C:\\Users\\law5_\\Desktop\\MDPerTool_GUI\\mdpertool\\output\\dis_protein_velocities.xml'
steps = range(0, 10, 1)  # Example steps to process

significant_atoms = track_significant_atoms(ref_file_path, dis_file_path, steps)
print(significant_atoms)

"""

"""
import xml.etree.ElementTree as ET
import numpy as np


# Example usage:
ref_file_path = 'C:\\Users\\law5_\\Desktop\\MDPerTool_GUI\\mdpertool\\output\\ref_protein_velocities1.xml'
pert_file_path = 'C:\\Users\\law5_\\Desktop\\MDPerTool_GUI\\mdpertool\\output\\dis_protein_velocities1.xml'

def parse_velocity_file(filename):
    tree = ET.parse(filename)
    root = tree.getroot()
    velocities = {}

    for velocity in root.findall('Velocity'):
        step = int(velocity.get('step', 0))
        atom_velocities = []
        for atom in velocity.findall('Atom'):
            x = float(atom.get('x', 0))
            y = float(atom.get('y', 0))
            z = float(atom.get('z', 0))
            atom_velocities.append([x, y, z])
        velocities[step] = np.array(atom_velocities)

    return velocities

def calculate_velocity_differences(ref_velocities, pert_velocities):
    steps = set(ref_velocities.keys()).union(set(pert_velocities.keys()))
    velocity_differences = {}

    for step in steps:
        ref_data = ref_velocities.get(step, np.zeros_like(next(iter(ref_velocities.values()))))
        pert_data = pert_velocities.get(step, np.zeros_like(next(iter(ref_velocities.values()))))
        velocity_differences[step] = pert_data - ref_data

    return velocity_differences

def count_affected_atoms(velocity_differences, threshold):
    affected_counts = {}
    for step, diffs in velocity_differences.items():
        magnitudes = np.linalg.norm(diffs, axis=1)
        affected_atoms = np.sum(magnitudes > threshold)
        affected_counts[step] = affected_atoms
    return affected_counts

threshold = 0.01  # Set your threshold value

# Read the data
ref_velocities = parse_velocity_file(ref_file_path)
pert_velocities = parse_velocity_file(pert_file_path)

# Calculate velocity differences
velocity_differences = calculate_velocity_differences(ref_velocities, pert_velocities)

# Count affected atoms
affected_atom_counts = count_affected_atoms(velocity_differences, threshold)

# Print results
for step, count in affected_atom_counts.items():
    print(f"Step {step}: Number of affected atoms: {count}")

"""

import os
import io
import tempfile
from pdbfixer import PDBFixer
from openmm.app import PDBFile


def fetched_pdb_fix(file_pathway, output_path=None, ph=7, chains_to_remove=None):
    """
    Args:
        :param file_pathway: Pathway for manipulating your fetched PDB files
        :param chains_to_remove: Selected chains will be deleted
        :param ph: Selected pH value will be applied to the structure's Hydrogens
    Returns:
        :param output_path: The manipulated PDB file will return as full path if specified
                            otherwise will return the already existing path
    """

    def remove_hetatoms_from_pdb(file_content):
        """
        Remove all HETATOM records from the PDB file content.

        Args:
            :param file_content: Content of the input PDB file as a string
        Returns:
            :return: Content of the PDB file with HETATOMs removed as a string
        """
        lines = file_content.splitlines()
        filtered_lines = [line for line in lines if not line.startswith("HETATM")]
        return '\n'.join(filtered_lines)

    # Read the original PDB file
    with open(file_pathway, 'r') as infile:
        pdb_content = infile.read()

    # Create a temporary file for the original PDB content
    with tempfile.NamedTemporaryFile(delete=False, mode='w', suffix='.pdb') as temp_file:
        temp_file.write(pdb_content)
        temp_file_path = temp_file.name

    print("Creating PDBFixer...")
    fixer = PDBFixer(temp_file_path)
    print("Finding missing residues...")

    if chains_to_remove is not None:
        print("toDelete: %s" % chains_to_remove)
        fixer.removeChains(chainIds=chains_to_remove)

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

    # Remove HETATOM records directly from the PDB content
    cleaned_pdb_content = remove_hetatoms_from_pdb(pdb_content)

    # Create a temporary file for the cleaned PDB content
    with tempfile.NamedTemporaryFile(delete=False, mode='w', suffix='.pdb') as temp_file:
        temp_file.write(cleaned_pdb_content)
        cleaned_pdb_path = temp_file.name

    # Read the cleaned PDB content with HETATOM removed
    fixer = PDBFixer(cleaned_pdb_path)

    print("Adding missing atoms...")
    fixer.findMissingResidues()
    print("Adding missing hydrogens...")
    fixer.addMissingHydrogens(pH=ph)

    print("Writing PDB file...")

    if output_path:
        path = os.path.join(output_path, f"{os.path.basename(file_pathway).split('.')[0]}_fixed_ph{ph}.pdb")
        with open(path, "w") as outfile:
            PDBFile.writeFile(fixer.topology, fixer.positions, outfile, keepIds=True)
        return path

    if not output_path:
        new_outpath = os.path.join(os.path.dirname(file_pathway),
                                   f"{os.path.basename(file_pathway).split('.')[0]}_fixed_ph{ph}.pdb")
        with open(new_outpath, "w") as outfile:
            PDBFile.writeFile(fixer.topology, fixer.positions, outfile, keepIds=True)
        return new_outpath


pdb = 'C:\\Users\\law5_\\Desktop\\MDPerTool_GUI\\mdpertool\\Download\\2j0w.pdb'
out = 'C:\\Users\\law5_\\Desktop\\MDPerTool_GUI\\mdpertool\\Download'
fetched_pdb_fix(pdb, output_path=out, ph=7, chains_to_remove=['A'])