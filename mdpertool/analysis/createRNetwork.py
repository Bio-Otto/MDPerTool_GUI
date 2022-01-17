import networkx as nx
from Bio.PDB.PDBParser import PDBParser
# from Bio.PDB import Residue
from PySide2.QtCore import QObject, Signal
from prody import *
# from collections import OrderedDict
import math
import matplotlib.pyplot as plt
# import csv
# import pandas
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
from analysis.multiproc_net_calc import Calc_Net_Worker
from .pdbsum_conservation_puller import *


def check_dissipated_residues_coordinates(pdb_file, querry_res_list, atom_name='CA'):
    dissipation_coordinates = {}
    p = PDBParser()
    structure = p.get_structure('prot', pdb_file)
    for model in structure:
        for chain in model:
            for residue in chain:
                ResId = residue.get_id()[1]
                ResName = residue.get_resname()
                Residue_label = ResName + str(ResId)

                if Residue_label in querry_res_list:
                    for atom in residue:
                        if atom.name == atom_name:
                            dissipation_coordinates[Residue_label] = atom.get_coord().tolist()

    return dissipation_coordinates


def check_shortest_path_residue_coordinates(pdb_file, querry_res_list, atom_name='CA'):
    shortest_path_res_coordinates = {}
    p = PDBParser()
    structure = p.get_structure('prot', pdb_file)
    for model in structure:
        for chain in model:
            for residue in chain:
                ResId = residue.get_id()[1]
                ResName = residue.get_resname()
                Residue_label = ResName + str(ResId)

                if Residue_label in querry_res_list:
                    for atom in residue:
                        if atom.name == atom_name:
                            shortest_path_res_coordinates[Residue_label] = atom.get_coord().tolist()

    return shortest_path_res_coordinates


def Pymol_Visualize_Path(graph, pdb_file):
    intersection_node_list = graph.nodes()  # GET NODES IN INTERSECTION GRAPH
    intersection_edge_list = graph.edges()  # GET EDGES IN INTERSECTION GRAPH
    coords = check_dissipated_residues_coordinates(pdb_file=pdb_file,
                                                   querry_res_list=intersection_node_list,
                                                   atom_name='CA')

    arrows_cordinates = []
    for edges in intersection_edge_list:
        arrows_cordinates.append((coords[edges[0]], coords[edges[1]]))

    return arrows_cordinates, intersection_node_list


def Shortest_Path_Visualize(pdb_file, selected_path):
    shortest_path_res_coordinates = check_shortest_path_residue_coordinates(pdb_file, querry_res_list=selected_path)

    arrows_cordinates = []
    for i in range(len(selected_path) - 1):
        arrows_cordinates.append(
            (shortest_path_res_coordinates[selected_path[i]], shortest_path_res_coordinates[selected_path[i + 1]]))

    return arrows_cordinates


def get_distance(array1, array2):
    x1, y1, z1 = float(array1[0]), float(array1[1]), float(array1[2])
    x2, y2, z2 = float(array2[0]), float(array2[1]), float(array2[2])
    return math.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2 + (z1 - z2) ** 2)


def get_residues(pdb):
    # gets the residue list from the pdb using PDBParser
    parser = PDBParser()
    structure = parser.get_structure('prot', pdb)
    # residue list in pdb
    residueList = list(structure.get_residues())
    # temp list for erasing the unnecessary atoms (water, hetatom etc.)
    resId_List = []
    garbageList = []
    for res in residueList:
        ResId = res.get_id()[1]
        ResName = res.get_resname()
        resId_List.append(ResName + str(ResId))

        if res.get_id()[0] != ' ':
            garbageList.append(res)
    # erase the residues with false id
    for res in garbageList:
        residueList.remove(res)

    return residueList, resId_List


def within_cutoff(res1, res2, distance_cutoff, just_CA):
    '''
    res1: BioPython PDB module Residue object.
    res2: BioPython PDB module Residue object.
    Calculates whether two residues have at least an atom pair within the specified cutoff distance.
    returns: True if residues are within cutoff distance, False if they are not.
    '''
    if just_CA:
        for atom1 in res1:
            if atom1.name == 'CA':
                atom1coords = atom1.get_coord()
                for atom2 in res2:
                    if atom2.name == 'CA':
                        atom2coords = atom2.get_coord()
                        distance = get_distance(atom1coords, atom2coords)
                        if distance <= distance_cutoff:
                            return True
    if not just_CA:
        for atom1 in res1:
            atom1coords = atom1.get_coord()
            for atom2 in res2:
                atom2coords = atom2.get_coord()
                distance = get_distance(atom1coords, atom2coords)
                if distance <= distance_cutoff:
                    return True
    return False


def createRNetwork(pdb, cutoff, reTimeFile, outputFileName, write_out, out_directory, CA_on=True):
    """
    residues: Biopython.PDB Residue objects representing the amino acids of a protein.
    distance_cutoff: Distance cutoff in Angstroms for edge (interaction) inclusion
    between two nodes (amino acids).
    """
    structure_res_list, res_list = get_residues(pdb=pdb)
    reTimeList = list()
    distance_cutoff = float(cutoff)
    network = nx.Graph()

    # ==> COLLECT RESPONSE TIMES FROM FILE
    with open(reTimeFile, 'r') as openReTimes:
        reTimes = openReTimes.readlines()
        for j in range(len(reTimes)):
            reTimeList.append(float(reTimes[j]))

    len_of_retimes_on_file = len(reTimes)

    # ==> ADD NODES AND ITS PROPERTIES TO GRAPH
    for res in range(len(structure_res_list)):
        node = structure_res_list[res]
        res_pos = node['CA'].get_coord()

        if not network.has_node(res_list[res]):
            network.add_node(res_list[res], posx=(str(res_pos[0])), posy=(str(res_pos[1])), posz=(str(res_pos[2])),
                             retime=(reTimeList[res]))
        else:
            network.add_node(res_list[res] + 'X', posx=(str(res_pos[0])), posy=(str(res_pos[1])),
                             posz=(str(res_pos[2])), retime=(reTimeList[res]))
            res_list[res] = res_list[res] + 'X'
            print("Residue of the same name was found while creating the RRI network. "
                  "Residue name changed from %s to %s." % (res_list[res][:-1], res_list[res]))

    # ==> ADD EDGES ACCORDING TO RESIDUE RESPONSE TIME AND CUT-OFF DISTANCE
    for i in range(len(structure_res_list)):
        residue1 = structure_res_list[i]
        for j in range(len(structure_res_list)):
            residue2 = structure_res_list[j]
            if res_list[i] != res_list[j]:
                if within_cutoff(residue1, residue2, distance_cutoff, just_CA=CA_on):
                    network.add_edge(res_list[i], res_list[j])

    if write_out:
        nx.write_gml(network, os.path.join(out_directory, outputFileName))

    return network, structure_res_list, res_list, len_of_retimes_on_file


def pairNetworks(network, source, target, pairNetworkName, write_out, out_directory, progress_callback,
                 node_threshold=None):
    targetRetime = network.nodes[target]['retime']
    sourceRetime = network.nodes[source]['retime']
    nodes = copy.deepcopy(network.nodes())
    log = ''

    # ==> FILTER OUT RESIDUES WITH A RESIDUE RESPONSE TIME GREATER THAN THE TARGET RESIDUE RESPONSE TIME
    for x in range(len(nodes)):
        nodeRetime = network.nodes[list(nodes)[x]]['retime']
        if nodeRetime > targetRetime:
            if list(nodes)[x] != target:
                network.remove_node(list(nodes)[x])

    # ==> COLLECT RESIDUE NAMES FROM Biopython.PDB
    res_list = []
    for i in network.nodes():
        res_list.append(str(i))

    clean_network = nx.DiGraph()
    # Adding Nodes
    for i in network.nodes():
        clean_network.add_node(i)

    # Adding Edges
    for i in range(len(res_list)):
        residue1 = res_list[i]
        retime1 = network.nodes[residue1]['retime']
        for j in range(len(res_list)):
            residue2 = res_list[j]
            if residue1 != residue2:
                retime2 = network.nodes[residue2]['retime']
                if network.has_edge(residue1, residue2):
                    if retime1 < retime2:
                        clean_network.add_edge(residue1, residue2)
                    if retime1 > retime2:
                        clean_network.add_edge(residue2, residue1)
                    if retime1 == retime2:
                        if (retime1 and retime2) != 0:
                            clean_network.add_edge(residue1, residue2)
                        if (retime1 and retime2) == 0:
                            if not clean_network.has_edge(residue1, residue2):
                                if residue1 == source:
                                    clean_network.add_edge(residue1, residue2)
                                if residue2 == source:
                                    clean_network.add_edge(residue2, residue1)

    # ==> FILTER OUT RESIDUES WITH IN-DEGREE EQUAL TO 0, EXCEPT SOURCE RESIDUE
    while clean_network.has_node(target):
        res_list = [str(i) for i in clean_network.nodes()]
        total_node_number = len(clean_network.nodes())
        junknodelist = list()

        for i in range(len(res_list)):
            chunk = res_list[i]
            if clean_network.in_degree(chunk) == 0:
                if chunk != source:  # and chunk != target
                    junknodelist.append(chunk)

        for junknode in junknodelist:
            clean_network.remove_node(junknode)

        total_node_number_2 = len(clean_network.nodes())
        if total_node_number == total_node_number_2:
            print('End of in-degree check')
            break

    if clean_network.has_node(target):
        # ==> FILTER OUT RESIDUES WITH OUT-DEGREE EQUAL TO 0, EXCEPT TARGET RESIDUE
        while True:
            res_list = [str(i) for i in clean_network.nodes()]
            total_node_number = len(clean_network.nodes())
            junknodelist2 = list()

            for i in range(len(res_list)):
                chunk = res_list[i]
                if clean_network.out_degree(chunk) == 0:
                    if chunk != target:  # and chunk != source
                        junknodelist2.append(chunk)

            for junknode in junknodelist2:
                clean_network.remove_node(junknode)

            total_node_number_2 = len(clean_network.nodes())
            if total_node_number == total_node_number_2:
                print('End of out-degree check')
                break

        nx.write_gml(clean_network, '%s' % pairNetworkName)

    else:
        print("YOUR TARGET RESIDUE ALREADY DELETED")

    if node_threshold is None and len(clean_network.nodes()) > 0:
        if clean_network.has_node(target):
            log = 'Source: %s  Target: %s\nTotal node number of source-target pair network is : %s' % (
                source, target, len(clean_network.nodes()))
            if write_out:
                nx.write_gml(clean_network, os.path.join(out_directory, pairNetworkName))

    if (node_threshold is not None) and (len(clean_network.nodes()) > node_threshold):
        if clean_network.has_node(target):
            log = 'Source: %s  Target: %s\nTotal node number of source-target pair network is : %s' % (
                source, target, len(clean_network.nodes()))
            if write_out:
                nx.write_gml(clean_network, os.path.join(out_directory, pairNetworkName))

    progress_callback.emit([clean_network, log])
    return clean_network, str(log)


class Multi_Task_Engine(object):

    def __init__(self, pdb_file, cutoff, reTimeFile, source, write_outputs, output_directory, node_threshold=None,
                 CA_on=True, outputFileName='111.gml', conserv_thresh=0.0, pdb_id=None):

        self.Work = []
        self.pdb_file = pdb_file
        self.cutoff = cutoff
        self.reTimeFile = reTimeFile
        self.outputFileName = outputFileName
        self.source = source
        self.write_outputs = write_outputs
        self.output_directory = output_directory
        self.node_threshold = node_threshold
        self.network = None
        self.resId_List = None
        self.CA_on = CA_on
        self.conserv_thresh = conserv_thresh
        self.pdb_id = pdb_id

    def calculate_general_network(self):
        try:
            self.network, residue_list, self.resId_List, len_of_reTimes = createRNetwork(pdb=self.pdb_file,
                                                                                         cutoff=self.cutoff,
                                                                                         reTimeFile=self.reTimeFile,
                                                                                         CA_on=self.CA_on,
                                                                                         outputFileName=self.outputFileName,
                                                                                         write_out=self.write_outputs,
                                                                                         out_directory=self.output_directory)
            return self.network, self.resId_List, len_of_reTimes

        except Exception as Error:
            print(Error)

    def run_pairNet_calc(self, target):
        print("====================================")
        print("TARGETS: ", target)
        print("====================================")
        try:

            for i in target:
                self.Work.append(
                    Calc_Net_Worker(pairNetworks, network=copy.deepcopy(self.network), source=self.source,
                                    target=i,
                                    pairNetworkName='%s_%s.gml' % (self.source, i),
                                    node_threshold=self.node_threshold,
                                    write_out=self.write_outputs,
                                    out_directory=self.output_directory))

        except IndexError as Err:
            print(Err)


