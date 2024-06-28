import networkx as nx
from Bio.PDB.PDBParser import PDBParser
from PySide2.QtCore import QObject, Signal
from prody import *
# from collections import OrderedDict
import math
import matplotlib.pyplot as plt
# import csv
import numpy as np
import argparse
from PySide2 import QtCore, QtWidgets
from src.PyMolWidget import PymolQtWidget
from .VisJS_Widget import VisJS_QtWidget
import multiprocessing as mp
from .pdbsum_conservation_puller import *
from Bio.PDB import Residue
from collections import OrderedDict
import math
import matplotlib.pyplot as plt
import csv
import pandas
import copy
import numpy as np
import os
# from analysis.multiproc_net_calc import Calc_Net_Worker
from analysis.multiproc_net_calc import CalcNetWorker
from .pdbsum_conservation_puller import *

import networkx as nx
import copy
import os
from typing import Callable, Optional, Tuple, List


def get_residue_coordinates(pdb_file, query_res_list, atom_name='CA'):
    """
    Get coordinates of specified atoms in the queried residues from a PDB file.

    Parameters:
    pdb_file (str): Path to the PDB file.
    query_res_list (list): List of residue labels to query (e.g., ['ALA12', 'GLY45']).
    atom_name (str): Name of the atom to get coordinates for (default is 'CA').

    Returns:
    dict: Dictionary with residue labels as keys and atom coordinates as values.
    """
    coordinates = {}
    p = PDBParser()
    structure = p.get_structure('prot', pdb_file)
    for model in structure:
        for chain in model:
            for residue in chain:
                res_id = residue.get_id()[1]
                res_name = residue.get_resname()
                residue_label = res_name + str(res_id)

                if residue_label in query_res_list:
                    try:
                        atom = residue[atom_name]
                        coordinates[residue_label] = atom.get_coord().tolist()
                    except KeyError:
                        print(f"Atom {atom_name} not found in residue {residue_label}")
    return coordinates


def check_dissipated_residues_coordinates(pdb_file, query_res_list, atom_name='CA'):
    """
    Get coordinates of specified atoms in the dissipated residues from a PDB file.

    Parameters:
    pdb_file (str): Path to the PDB file.
    query_res_list (list): List of dissipated residue labels to query.
    atom_name (str): Name of the atom to get coordinates for (default is 'CA').

    Returns:
    dict: Dictionary with residue labels as keys and atom coordinates as values.
    """
    return get_residue_coordinates(pdb_file, query_res_list, atom_name)


def check_shortest_path_residue_coordinates(pdb_file, query_res_list, atom_name='CA'):
    """
    Get coordinates of specified atoms in the shortest path residues from a PDB file.

    Parameters:
    pdb_file (str): Path to the PDB file.
    query_res_list (list): List of shortest path residue labels to query.
    atom_name (str): Name of the atom to get coordinates for (default is 'CA').

    Returns:
    dict: Dictionary with residue labels as keys and atom coordinates as values.
    """
    return get_residue_coordinates(pdb_file, query_res_list, atom_name)


def Pymol_Visualize_Path(graph, pdb_file):
    """
    Get coordinates for visualizing a path in PyMOL from a graph and a PDB file.

    Parameters:
    graph (nx.Graph): The graph representing the path.
    pdb_file (str): Path to the PDB file.

    Returns:
    tuple: A tuple containing:
        - List of arrow coordinates for PyMOL visualization.
        - List of nodes in the intersection graph.
    """
    intersection_node_list = list(graph.nodes)
    intersection_edge_list = list(graph.edges)
    coords = check_dissipated_residues_coordinates(pdb_file, intersection_node_list, 'CA')

    arrow_coordinates = [(coords[edge[0]], coords[edge[1]]) for edge in intersection_edge_list]

    return arrow_coordinates, intersection_node_list


def Shortest_Path_Visualize(pdb_file, selected_path):
    """
    Get coordinates for visualizing the shortest path in PyMOL from a PDB file.

    Parameters:
    pdb_file (str): Path to the PDB file.
    selected_path (list): List of residue labels representing the selected path.

    Returns:
    list: List of arrow coordinates for PyMOL visualization.
    """
    coords = check_shortest_path_residue_coordinates(pdb_file, selected_path, 'CA')

    arrow_coordinates = [(coords[selected_path[i]], coords[selected_path[i + 1]]) for i in
                         range(len(selected_path) - 1)]

    return arrow_coordinates


def get_distance(coord1, coord2):
    """
    Calculate the Euclidean distance between two coordinates.

    Parameters:
    coord1 (array-like): Coordinates of the first point.
    coord2 (array-like): Coordinates of the second point.

    Returns:
    float: The Euclidean distance between the two points.
    """
    return np.linalg.norm(coord1 - coord2)


def get_residues(pdb):
    """
    Get residue list from the pdb using PDBParser.

    Parameters:
    pdb (str): Path to the PDB file.

    Returns:
    tuple: A tuple containing:
        - list of residues (Bio.PDB.Residue.Residue)
        - list of residue identifiers (str)
    """
    parser = PDBParser()
    structure = parser.get_structure('prot', pdb)
    residue_list = [res for res in structure.get_residues() if res.get_id()[0] == ' ']
    res_id_list = [res.get_resname() + str(res.get_id()[1]) for res in residue_list]
    return residue_list, res_id_list


def within_cutoff(res1, res2, distance_cutoff, just_CA):
    """
    Calculate whether two residues have at least an atom pair within the specified cutoff distance.

    Parameters:
    res1 (Bio.PDB.Residue.Residue): The first residue.
    res2 (Bio.PDB.Residue.Residue): The second residue.
    distance_cutoff (float): The cutoff distance in Angstroms.
    just_CA (bool): Whether to consider only CA atoms.

    Returns:
    bool: True if the residues are within the cutoff distance, False otherwise.
    """
    if just_CA:
        ca1 = res1['CA'].get_coord()
        ca2 = res2['CA'].get_coord()
        return get_distance(ca1, ca2) <= distance_cutoff

    for atom1 in res1:
        for atom2 in res2:
            if get_distance(atom1.get_coord(), atom2.get_coord()) <= distance_cutoff:
                return True

    return False


def load_response_times(reTimeFile):
    """
    Load response times from the specified file.

    Parameters:
    reTimeFile (str): Path to the file containing response times.

    Returns:
    list: A list of response times (float).
    """
    reTimeList = []
    with open(reTimeFile, 'r') as file:
        reTimeList = [float(line.strip()) for line in file]
    return reTimeList


def createRNetwork(pdb, cutoff, reTimeFile, outputFileName, write_out, out_directory, CA_on=True):
    """
    Create a residue interaction network based on a distance cutoff and residue response times.

    Parameters:
    pdb (str): Path to the PDB file.
    cutoff (float): Distance cutoff in Angstroms for edge inclusion between residues.
    reTimeFile (str): Path to the file containing residue response times.
    outputFileName (str): Name of the output file to save the network.
    write_out (bool): Whether to write the output network to a file.
    out_directory (str): Directory to write the output file to.
    CA_on (bool): Whether to consider only CA atoms for distance calculations.

    Returns:
    tuple: A tuple containing:
        - network (nx.Graph): The residue interaction network.
        - structure_res_list (list): List of Bio.PDB.Residue.Residue objects.
        - res_list (list): List of residue identifiers.
        - len_of_retimes_on_file (int): The number of response times loaded from the file.
    """
    try:
        structure_res_list, res_list = get_residues(pdb)
        reTimeList = load_response_times(reTimeFile)
        distance_cutoff = float(cutoff)
        network = nx.Graph()

        len_of_retimes_on_file = len(reTimeList)

        # Add nodes and properties to the graph
        for idx, res in enumerate(structure_res_list):
            res_pos = res['CA'].get_coord()
            node_name = res_list[idx]

            if network.has_node(node_name):
                node_name += 'X'
                res_list[idx] = node_name
                print(f"Residue of the same name was found. Residue name changed to {node_name}.")

            network.add_node(node_name, posx=str(res_pos[0]), posy=str(res_pos[1]), posz=str(res_pos[2]),
                             retime=reTimeList[idx])

        # Add edges according to residue response time and cutoff distance
        for i, res1 in enumerate(structure_res_list):
            for j, res2 in enumerate(structure_res_list):
                if i != j and within_cutoff(res1, res2, distance_cutoff, just_CA=CA_on):
                    network.add_edge(res_list[i], res_list[j])

        if write_out:
            nx.write_gml(network, os.path.join(out_directory, outputFileName))

        return network, structure_res_list, res_list, len_of_retimes_on_file

    except Exception as e:
        print(f"An error occurred: {e}")
        return None, None, None, 0


def filter_pair_network(network, source, target, node_threshold=None):
    """
    Filters and creates a pair network between a source and target residue.

    Args:
        network (networkx.DiGraph): The input network.
        source (str): The source residue.
        target (str): The target residue.
        node_threshold (int, optional): The maximum number of nodes allowed in the pair network.

    Returns:
        networkx.DiGraph: The filtered pair network.
    """
    # Filter out residues with a residue response time greater than the target residue response time
    target_retime = network.nodes[target]['retime']
    clean_network = nx.DiGraph()
    for node, data in network.nodes(data=True):
        if data['retime'] <= target_retime or node == target:
            clean_network.add_node(node, **data)

    # Add edges to the clean network
    for u, v, data in network.edges(data=True):
        u_retime = network.nodes[u]['retime']
        v_retime = network.nodes[v]['retime']
        if u_retime < v_retime:
            clean_network.add_edge(u, v, **data)
        elif u_retime > v_retime:
            clean_network.add_edge(v, u, **data)
        elif (u_retime and v_retime) != 0:
            clean_network.add_edge(u, v, **data)
        elif (u_retime and v_retime) == 0 and (u == source or v == source):
            clean_network.add_edge(u, v, **data)

    # Filter out residues with in-degree equal to 0, except the source residue
    while True:
        junk_nodes = [node for node in clean_network.nodes() if clean_network.in_degree(node) == 0 and node != source]
        if not junk_nodes:
            break
        clean_network.remove_nodes_from(junk_nodes)

    # Filter out residues with out-degree equal to 0, except the target residue
    while True:
        junk_nodes = [node for node in clean_network.nodes() if clean_network.out_degree(node) == 0 and node != target]
        if not junk_nodes:
            break
        clean_network.remove_nodes_from(junk_nodes)

    # Apply the node threshold if specified
    if node_threshold is not None and len(clean_network.nodes()) > node_threshold:
        return None
    else:
        return clean_network


def pairNetworks(network, source, target, pairNetworkName, write_out, out_directory, progress_callback, node_threshold=None):
    """
    Filters and creates a pair network between a source and target residue.

    Args:
        network (networkx.DiGraph): The input network.
        source (str): The source residue.
        target (str): The target residue.
        pairNetworkName (str): The name of the pair network.
        write_out (bool): Whether to write the pair network to a file.
        out_directory (str): The directory to write the pair network file.
        progress_callback (callable): A callback function to report progress.
        node_threshold (int, optional): The maximum number of nodes allowed in the pair network.

    Returns:
        networkx.DiGraph: The filtered pair network.
        str: A log message with information about the pair network.
    """
    clean_network = filter_pair_network(network, source, target, node_threshold)
    if clean_network is not None:
        log = f'Source: {source}  Target: {target}\nTotal node number of source-target pair network is: {len(clean_network.nodes())}'
        if write_out:
            nx.write_gml(clean_network, os.path.join(out_directory, pairNetworkName))
    else:
        log = f'Source: {source}  Target: {target}\nPair network not created due to node threshold.'

    progress_callback.emit([clean_network, log])
    return clean_network, log


class MultiTaskEngine:
    """
    A class to handle multiple tasks related to network calculations based on PDB files.
    """

    def __init__(self, pdb_file, cutoff, re_time_file, source, write_outputs, output_directory, node_threshold=None,
                 ca_on=True, output_file_name='111.gml', conserv_thresh=0.0, pdb_id=None):
        """
        Initialize the MultiTaskEngine with the necessary parameters.

        :param pdb_file: Path to the PDB file
        :param cutoff: Cutoff value for network calculations
        :param re_time_file: File for re-time calculations
        :param source: Source for network calculations
        :param write_outputs: Boolean to decide whether to write outputs
        :param output_directory: Directory to write outputs
        :param node_threshold: Threshold for nodes in the network
        :param ca_on: Boolean to decide if CA atoms should be considered
        :param output_file_name: Name of the output file
        :param conserv_thresh: Conservation threshold
        :param pdb_id: PDB ID
        """
        self.work = []
        self.pdb_file = pdb_file
        self.cutoff = cutoff
        self.re_time_file = re_time_file
        self.output_file_name = output_file_name
        self.source = source
        self.write_outputs = write_outputs
        self.output_directory = output_directory
        self.node_threshold = node_threshold
        self.network = None
        self.res_id_list = None
        self.ca_on = ca_on
        self.conserv_thresh = conserv_thresh
        self.pdb_id = pdb_id

    def calculate_general_network(self):
        """
        Calculate the general network based on the provided PDB file and parameters.

        :return: A tuple containing the network, residue ID list, and length of re-times.
        """
        try:
            self.network, residue_list, self.res_id_list, len_of_re_times = createRNetwork(
                pdb=self.pdb_file,
                cutoff=self.cutoff,
                reTimeFile=self.re_time_file,
                CA_on=self.ca_on,
                outputFileName=self.output_file_name,
                write_out=self.write_outputs,
                out_directory=self.output_directory
            )
            return self.network, self.res_id_list, len_of_re_times

        except Exception as error:
            print(f"Error calculating general network: {error}")

    def run_pair_network_calculation(self, targets):
        """
        Run pair network calculations for the given targets.

        :param targets: List of target nodes for network calculations
        """
        # print("====================================")
        # print("TARGETS:", targets)
        # print("====================================")
        try:
            for target in targets:
                # print(target)
                self.work.append(
                    CalcNetWorker(
                        pairNetworks,
                        network=copy.deepcopy(self.network),
                        source=self.source,
                        target=target,
                        pairNetworkName=f'{self.source}_{target}.gml',
                        node_threshold=self.node_threshold,
                        write_out=self.write_outputs,
                        out_directory=self.output_directory
                    )
                )

        except IndexError as error:
            print(f"Index error during pair network calculation: {error}")
