import os
import time
import sys
import csv
import copy
from pathlib import Path
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
from PyMolWidget import PymolQtWidget
from .VisJS_Widget import VisJS_QtWidget
import multiprocessing as mp
from .pdbsum_conservation_puller import *
from analysis.multiproc_net_calc import Calc_Net_Worker


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


def createRNetwork(pdb, cutoff, reTimeFile, outputFileName, write_out, out_directory, verbose=False):
    '''creates a residue network having edges within given cutoff range'''

    global posx, posy, posz

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

    def calculate_distance(coords1, coords2):
        # calculates the distance between two atoms
        xDistance = ((coords1[0]) - (coords2[0])) ** 2
        yDistance = ((coords1[1]) - (coords2[1])) ** 2
        zDistance = ((coords1[2]) - (coords2[2])) ** 2
        xyzDistance = math.sqrt((xDistance) + (yDistance) + (zDistance))

        return xyzDistance

    # def withinCutoff(res1, res2, cutoff, verbose=False):
    # 	#decides whether the atom is within the cutoff range or not. If it is, returns True.
    # 	for atom1 in res1:
    # 		atom1coords = atom1.get_coord()
    # 		for atom2 in res2:
    # 			atom2coords = atom2.get_coord()
    # 			distance = calculate_distance(atom1coords, atom2coords)
    # 			if (distance <= cutoff):
    # 				return True

    # print('Cutoff is %d' % (cutoff))

    # if verbose:
    #     print('Initializing...')

    # start creating network
    residueList, resId_List = get_residues(pdb)
    network = nx.DiGraph()

    # create a nodelist for proper naming
    nodeList = []

    for chunk in range(0, len(residueList)):
        node = residueList[chunk]
        nodeNumber = str(node.id[1])
        nodeID = node.get_resname() + nodeNumber
        nodeList.append(nodeID)
    # extract the positions of calpas and add nodes with position and retime attributes.
    openReTimes = open(reTimeFile, 'r')
    reTimeList = list()
    reTimes = openReTimes.read().splitlines()

    for j in range(0, len(reTimes)):
        chunk = reTimes[j]
        reTimeList.append(chunk)

    for res in range(0, len(residueList)):
        node = residueList[res]

        for atom in node:
            if atom.name == 'CA':
                posx, posy, posz = atom.get_coord()

        network.add_node(nodeList[res], posx=('%f' % posx), posy=('%f' % posy), posz=('%f' % posz),
                         retime=(reTimeList[res]))

    # if verbose:
    #     print('Adding nodes to network...')

    # print('Calculating distances...')
    # print('Assigning directions...')
    # add edges between nodes if they are within cutoff range and decide the direction using retimes of nodes
    # if the retime of node A is less than the node B direction is from A to B.
    for i in range(0, len(residueList)):
        residue1 = residueList[i]
        # coords1 = residue1['CA'].get_coord()

        for atom1 in residue1:
            if atom1.name == 'CA':
                coords1 = atom1.get_coord()

        retime1 = network.nodes[nodeList[i]]['retime']
        for j in range(i + 1, len(residueList)):
            residue2 = residueList[j]

            for atom2 in residue2:
                if atom2.name == 'CA':
                    coords2 = atom2.get_coord()

            retime2 = network.nodes[nodeList[j]]['retime']
            distance = calculate_distance(coords1, coords2)
            if distance <= cutoff:
                if retime1 <= retime2:
                    network.add_edge(nodeList[i], nodeList[j])
                else:
                    network.add_edge(nodeList[j], nodeList[i])

    # if verbose:
    #     print('Adding edges between nodes...')

    if write_out:
        nx.write_gml(network, os.path.join(out_directory, outputFileName))

    # if verbose:
    #     print('Network construction is complete.')
    # print network.node
    # position = nx.get_node_attributes(network, "pos")
    # nx.draw_networkx_nodes(network, position)
    # nx.draw_networkx_edges(network, position)
    # plt.show()
    return network, residueList, resId_List


def pairNetworks(network, source, target, pairNetworkName, write_out, out_directory, progress_callback,
                 node_threshold=None, verbose=True):
    sourceRetime = network.nodes['%s' % (source)]['retime']
    targetRetime = network.nodes['%s' % (target)]['retime']

    log = ''
    # newNetwork = nx.DiGraph()
    # print('Initializing source-target pair network...')

    nodes = list(network.nodes().keys())

    # print('Checking retime of nodes')

    for x in nodes:
        nodeRetime = network.nodes[x]['retime']
        if x != target:
            if nodeRetime > targetRetime:
                network.remove_node(x)
                # print('Removing node %s' % x)

    newNetwork = copy.deepcopy(network)

    while True:
        if target not in newNetwork.nodes().keys():
            pass
            # print("Target Residue Deleted In Degree Process !!!!!!!")

        total_node_number = len(newNetwork.nodes())

        for i in list(network.nodes().keys()):
            inDegree = newNetwork.in_degree(i)

            # print('Checking in-degree of %s' % i)
            if inDegree == 0:

                if i != source:
                    try:
                        # print('Checking in-degree of %s' % i)
                        newNetwork.remove_node(i)

                    except ValueError:
                        pass
                        # print("Indicated node (%s) already deleted" % i)

                # print('Indegree of %s is %s' % (i, inDegree))

        if total_node_number == len(newNetwork.nodes()):
            # print('End of in-degree check')
            break

    # print("%s" % len(newNetwork.nodes()))

    while True:
        total_node_number = len(newNetwork.nodes())

        for j in list(network.nodes().keys()):
            outDegree = newNetwork.out_degree(j)

            # print('Checking out-degree of %s' % j)
            if outDegree == 0:
                if j != target:
                    try:
                        newNetwork.remove_node(j)
                    except:
                        pass
                        # print("Indicated node (%s) already deleted" % j)

                # print('Outdegree of %s is %s' % (j, outDegree))

        total_node_number_2 = len(newNetwork.nodes())

        if total_node_number == total_node_number_2:
            # print('End of out-degree check')
            break

    if node_threshold is None and len(newNetwork.nodes()) > 0:
        if verbose:
            # print("Source: %s  Target: %s" % (source, target))
            # print('Total node number of source-target pair network is : ', len(newNetwork.nodes()))
            log = 'Source: %s  Target: %s\nTotal node number of source-target pair network is : %s' % (
                source, target, len(newNetwork.nodes()))
        if write_out:
            nx.write_gml(newNetwork, os.path.join(out_directory, pairNetworkName))

    if (node_threshold is not None) and (len(newNetwork.nodes()) > node_threshold):
        if verbose:
            # print("Source: %s  Target: %s" % (source, target))
            # print('Total node number of source-target pair network is : ', len(newNetwork.nodes()))
            log = 'Source: %s  Target: %s\nTotal node number of source-target pair network is : %s' % (
                source, target, len(newNetwork.nodes()))
        if write_out:
            nx.write_gml(newNetwork, os.path.join(out_directory, pairNetworkName))

    # print('Done! Ready to Visualize')
    progress_callback.emit([newNetwork, log])
    return newNetwork, str(log)


def getResNetwork(resContMap, CorrFile, pdb, outputFileName, verbose=True):
    def get_residues(pdb):

        # gets the residue list from the pdb using PDBParser

        parser = PDBParser()
        structure = parser.get_structure('prot', pdb)
        # residue list in pdb
        residueList = list(structure.get_residues())
        # temp list for erasing the unnecessary atoms (water, hetatom etc.)
        garbageList = []
        for res in residueList:
            if res.get_id()[0] != ' ':
                garbageList.append(res)
        # erase the residues with false id
        for res in garbageList:
            residueList.remove(res)
        return residueList

    residueList = get_residues(pdb)
    network = nx.Graph()
    nodeList = []
    sys = parsePDB(pdb)
    calphas = sys.select('calpha')
    resnums = calphas.getResnums().tolist()

    # create a nodelist for proper naming
    for chunk in range(0, len(residueList)):
        node = residueList[chunk]
        nodeNumber = str(node.id[1])
        nodeID = node.get_resname() + nodeNumber
        nodeList.append(nodeID)
        network.add_node(nodeList[chunk])

        print(network.nodes())
    # Load the CorrFile and determine the maximum edge weight
    resCorrMat = np.loadtxt(CorrFile)
    resCorrArray = np.squeeze(resCorrMat)
    maxResCorr = 0.4

    if verbose:
        print('Adding nodes to network...')

    # read Contact map and assign edges

    f = open(resContMap, 'r')
    lines = f.readlines()

    for i in range(1, len(lines)):
        line = lines[i].split()
        print(line)
        res1name = line[0].split(":")
        res1 = res1name[3] + res1name[1]
        res2name = line[2].split(":")
        res2 = res2name[3] + res2name[1]
        res1num = resnums.index(int(res1name[1]))
        res2num = resnums.index(int(res2name[1]))
        distance = float(resCorrMat[res1num, res2num])
        if distance >= maxResCorr:
            network.add_edge(res1, res2, {'distance': distance})
        else:
            continue

    nx.write_gml(network, '%s.gml' % outputFileName)
    if verbose:
        print('Network construction is complete.')

    return network


def analyzeNetwork(resContMap, CorrFile, pdb, outputFileName, source, verbose=False):
    N = getResNetwork(resContMap, CorrFile, pdb, outputFileName, True)

    if verbose:
        print('Calculating shortest paths from %s' % source)
        print('Calculating closeness centrality...')
        print('Calculating betweenness centrality...')

    shortest_paths = nx.single_source_shortest_path(N, source)
    close = nx.closeness_centrality(N)
    between = nx.betweenness_centrality(N)
    columns = list()

    for n in N:
        data = '%s, %3f, %3f' % (n, close[n], between[n])
        columns.append(data)
    if verbose:
        print('Writing files...')

    outputFile = open(outputFileName + 'analysis.csv', 'w')
    outputFile.write('Residue ID, Closeness Centrality,Betweenness Centrality\n')
    for i in range(0, len(columns)):
        outputFile.write(columns[i])
        outputFile.write('\n')
    outputFile.close()

    outputFile_2 = open(outputFileName + '_shortestpath.csv', 'w')
    for m in shortest_paths:
        path = shortest_paths[m]
        outputFile_2.write(m + ':')
        for k in range(0, len(path)):
            outputFile_2.write(path[k] + ',')
        outputFile_2.write('\n')
    outputFile_2.close()
    return N


def runTest(pdb, cutoff, reTimeFile, outputFileName, source, target, pairNetworkName):
    network, residue_list = createRNetwork(pdb, cutoff, reTimeFile, outputFileName, verbose=False, write_out=True,
                                           out_directory=os.getcwd())
    pairNetworks(network, source, target, pairNetworkName)
    return network


def main():
    parser = argparse.ArgumentParser(description="Run Test Aalyzing Module")
    parser.add_argument("-pdb", help="pdb input file", dest="pdb", type=str, required=True)
    parser.add_argument("-cutoff", help="cutoff input value", dest="cutoff", type=int, required=True)
    parser.add_argument("-reTimeFile", help="ReTimeFile input file", dest="reTimeFile", type=str, required=True)
    parser.add_argument("-out", help="output file", dest="outputFileName", type=str, required=True)
    parser.add_argument("-source", help="Source", dest="source", type=str, required=True)
    parser.add_argument("-target", help="Indicate Target", dest="target", type=str, required=True)
    parser.add_argument("-pairNetworkName", help="pair Network Name", dest="pairNetworkName", type=str, required=True)
    parser.set_defaults(func=runTest)
    args = parser.parse_args()
    args.func(args)


class Multi_Task_Engine(object):

    def __init__(self, pdb_file, cutoff, reTimeFile, source, write_outputs, output_directory, node_threshold=None,
                 verbose=True, outputFileName='111.gml', conserv_thresh=0.0, pdb_id=None):

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
        self.verbose = verbose
        self.conserv_thresh = conserv_thresh
        self.pdb_id = pdb_id

    def calculate_general_network(self):
        try:
            self.network, residue_list, self.resId_List = createRNetwork(pdb=self.pdb_file, cutoff=self.cutoff,
                                                                         reTimeFile=self.reTimeFile,
                                                                         verbose=self.verbose,
                                                                         outputFileName=self.outputFileName,
                                                                         write_out=self.write_outputs,
                                                                         out_directory=self.output_directory)
            return self.network, self.resId_List

        except Exception as Error:
            print(Error)

    def run_pairNet_calc(self, target):
        try:

            for i in target:

                self.Work.append(
                    Calc_Net_Worker(pairNetworks, network=copy.deepcopy(self.network), source=self.source,
                                    target=i,
                                    pairNetworkName='%s_%s.gml' % (self.source, target),
                                    node_threshold=self.node_threshold,
                                    verbose=self.verbose, write_out=self.write_outputs,
                                    out_directory=self.output_directory))
            # network, log = pairNetworks(network=self.calculate_general_network[0], source=self.source, target=target,
            #                             pairNetworkName='%s_%s.gml' % (self.source, target),
            #                             node_threshold=self.node_threshold,
            #                             verbose=self.verbose, write_out=self.write_outputs,
            #                             out_directory=self.output_directory)

            # return network, log
        except IndexError as Err:
            print(Err)


def intersection_of_directed_networks(graphs_list):
    len_of_nodes_on_list = [len(graph.nodes()) for graph in graphs_list]
    smallest_network_and_indices = min([(v, i) for i, v in enumerate(len_of_nodes_on_list)])

    R = copy.deepcopy(graphs_list[smallest_network_and_indices[1]])

    for cnt, graph_count in enumerate(graphs_list):
        if cnt != smallest_network_and_indices[1]:
            R.remove_nodes_from(n for n in graphs_list[smallest_network_and_indices[1]] if n not in graphs_list[cnt])
            R.remove_edges_from(
                e for e in graphs_list[smallest_network_and_indices[1]].edges if e not in graphs_list[cnt].edges)

    return R


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


def call_pymol_for_network_visualization(pdb_file, arrows_cordinates, intersection_node_list):
    if not QtWidgets.QApplication.instance():
        app = QtWidgets.QApplication(sys.argv)
    else:
        app = QtWidgets.QApplication.instance()

    window = PymolQtWidget()
    window.resize(1200, 800)
    window.loadMolFile(pdb_file)
    window.show_energy_dissipation(response_time_file_path=retime_file)

    for arrow_coord in arrows_cordinates:
        window.create_directed_arrows(atom1=arrow_coord[0], atom2=arrow_coord[1], radius=0.05,
                                      gap=0.4, hradius=0.4, hlength=0.8, color='green')

    for node in intersection_node_list:
        resID_of_node = int(''.join(list(filter(str.isdigit, node))))
        window.resi_label_add('resi ' + str(resID_of_node))

    # MAKE PYMOL VISUALIZATION BETTER
    window._pymol.cmd.set('cartoon_oval_length', 0.8)  # default is 1.20)
    window._pymol.cmd.set('cartoon_oval_width', 0.2)
    window._pymol.cmd.center(selection="all", state=0, origin=1, animate=0)
    window._pymol.cmd.zoom('all', buffer=0.0, state=0, complete=0)

    window.update()
    window.show()
    app.exec_()


def call_visJS_for_network_visualization(network):
    if not QtWidgets.QApplication.instance():
        app = QtWidgets.QApplication(sys.argv)
    else:
        app = QtWidgets.QApplication.instance()

    w = VisJS_QtWidget(network=network)
    return w, app

    # sys.exit(app.exec_())


"""
if __name__ == '__main__':  # PROTECT YOUR PROGRAM'S ENTRY POINT
    global intersection_graph, output_folder_directory

    # --------------------------------------------> INPUTS / START <-------------------------------------------------- #
    pdb = 'test/Coupled_Local_Motion/1UFP/2EB8.pdb'  # PDB file path --> "BOUND FORM OF STRUCTURE"
    cutoff = 10  # Add edges between nodes if they are within cutoff range
    retime_file = 'test/Coupled_Local_Motion/1UFP/P1/responseTimes_5.csv'  # Response time file path
    outputFileName = 'protein_general_cutoff_network.gml'  # Protein general graph according to cut off value
    source = 'SER92'  # One of the perturbed residues
    node_threshold = 20  # None or an Integer
    verbose_condition = True  # True or False
    target_residues = ['GLU54', 'GLY25', 'HIS64', 'THR67']  # ['PHE18', 'ASP83', 'TYR35', 'ALA50', 'ILE61'] # None or residue list

    use_conservation = True
    pdb_id = '2EB8'  # FREE OR BOUND FORM OF PDB
    chain = 'A'  # FOR PULLING CONSERVATION SCORES INDICATE CHAIN ID OF PDB
    conservation_threshold = 9
    save_conservation_scores = True

    visualize_on_PyMol = True  # Networkx Graph Visulization on PyMol
    visualize_on_VisJS = True  # Networkx Graph Visulization on VisJS
    create_output = True  # Supports True or False Conditions for creation of all networks (*.gml) on a folder
    just_visualize = False  # If you have already calculated network you can directly visualize it.

    # ---------------------------------------------- INPUTS / END <---------------------------------------------- #

    general_output_folder = os.path.join(os.getcwd(), 'OUTPUTS')
    Path(general_output_folder).mkdir(parents=True, exist_ok=True)

    if not just_visualize:
        mp.freeze_support()
        start_Time = time.time()  # GET START TIME
        pool = mp.Pool(mp.cpu_count())  # CREATE A POOL WITH MAX CPU THAT YOU HAVE IN YOUR MACHINE

    if create_output:
        folder_name = "output_%s" % source
        output_folder_directory = os.path.join(general_output_folder, folder_name)
        Path(output_folder_directory).mkdir(parents=True, exist_ok=True)

        if just_visualize:
            if os.path.exists(os.path.join(output_folder_directory, 'intersection_graph.gml')):
                print(os.path.join(output_folder_directory, 'intersection_graph.gml'))
                intersection_graph = nx.read_gml(os.path.join(output_folder_directory, 'intersection_graph.gml'))

                # ----------------------> CALL PYMOL FOR VISUALIZATION / START <---------------------- #
                if visualize_on_PyMol:
                    print(" PyMol Visualization")
                    arrows_cordinates, intersection_node_list = Pymol_Visualize_Path(graph=intersection_graph,
                                                                                     pdb_file=pdb)
                    call_pymol_for_network_visualization(pdb_file=pdb, arrows_cordinates=arrows_cordinates,
                                                         intersection_node_list=intersection_node_list)

                if visualize_on_VisJS:
                    visjs_engine, app = call_visJS_for_network_visualization(network=intersection_graph)
                    visjs_engine()
                    visjs_engine.show()
                    sys.exit(app.exec_())

                if visualize_on_PyMol is False and visualize_on_VisJS is False:
                    print("There is no active visualization choice. "
                          "Make one of the visualization program condition True")
                sys.exit()
                # -----------------------> CALL PYMOL FOR VISUALIZATION / END <----------------------- #
            else:
                print("There is no *.gml file for visualization.")
                sys.exit()

    try:
        if create_output:
            engine = Multi_Task_Engine(pdb_file=pdb, cutoff=cutoff, reTimeFile=retime_file, source=source,
                                       node_threshold=node_threshold, verbose=verbose_condition,
                                       outputFileName=outputFileName, write_outputs=create_output,
                                       output_directory=output_folder_directory)

        if not create_output:
            engine = Multi_Task_Engine(pdb_file=pdb, cutoff=cutoff, reTimeFile=retime_file, source=source,
                                       node_threshold=node_threshold, verbose=verbose_condition,
                                       outputFileName=outputFileName, write_outputs=create_output,
                                       output_directory=None)

        if target_residues is None:
            if use_conservation:
                res_IDs, con_scores = get_conservation_scores(pdb_id=pdb_id, chain_id=chain,
                                                              cutoff=conservation_threshold, bound_pdb=pdb)

                if save_conservation_scores:
                    rows = zip(res_IDs, con_scores)
                    with open(os.path.join(output_folder_directory, 'conservation_%s.csv' % pdb_id), "w",
                              newline='') as f:
                        writer = csv.writer(f)
                        for row in rows:
                            writer.writerow(row)

                net = pool.map(engine, res_IDs)

            if not use_conservation:
                resId_list = engine.calculate_general_network[1]
                net = pool.map(engine, resId_list)

        if target_residues is not None:
            if use_conservation:
                res_IDs, con_scores = get_conservation_scores(pdb_id=pdb_id, chain_id=chain,
                                                              cutoff=conservation_threshold, bound_pdb=pdb)
                if save_conservation_scores:
                    rows = zip(res_IDs, con_scores)
                    with open(os.path.join(output_folder_directory, 'conservation_%s.csv' % pdb_id), "w",
                              newline='') as f:
                        writer = csv.writer(f)
                        for row in rows:
                            writer.writerow(row)
                intersection_resIDs = set.intersection(set(res_IDs), set(target_residues))
                net = pool.map(engine, list(intersection_resIDs))

            if not use_conservation:
                net = pool.map(engine, target_residues)

        clean_graph_list = []
        if node_threshold is not None:
            for i in net:
                if len(i.nodes()) > node_threshold:
                    clean_graph_list.append(i)

        if node_threshold is None:
            for i in net:
                clean_graph_list.append(i)

        # CREATE AN INTERSECTION GRAPH AND WRITE TO GML FILE
        if len(clean_graph_list) > 0:
            intersection_graph = intersection_of_directed_networks(clean_graph_list)
            if create_output:
                nx.write_gml(intersection_graph, os.path.join(output_folder_directory, 'intersection_graph.gml'))

        else:
            print("There is no suitable Graph for your search parameters")

        # VISUALIZE DIRECTED NETWORK ON PYMOL
        if visualize_on_PyMol:
            try:
                arrows_cordinates, intersection_node_list = Pymol_Visualize_Path(graph=intersection_graph, pdb_file=pdb)

                # ----------------------> CALL PYMOL FOR VISUALIZATION / START <---------------------- #
                call_pymol_for_network_visualization(pdb_file=pdb, arrows_cordinates=arrows_cordinates,
                                                     intersection_node_list=intersection_node_list)
                # -----------------------> CALL PYMOL FOR VISUALIZATION / END <----------------------- #

            except Exception as error:
                print(error)

        if visualize_on_VisJS:
            try:

                # ----------------------> CALL VisJS FOR VISUALIZATION / START <---------------------- #
                visjs_engine, app = call_visJS_for_network_visualization(network=intersection_graph)
                visjs_engine()
                visjs_engine.show()
                sys.exit(app.exec_())
                # -----------------------> CALL VisJS FOR VISUALIZATION / END <----------------------- #

            except Exception as error:
                print(error)

    except Exception as e:
        exc_type, exc_obj, exc_tb = sys.exc_info()
        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        print(exc_type, fname, exc_tb.tb_lineno)

    finally:
        # CLOSE THE ALREADY OPENED POOL
        pool.close()
        pool.join()
        print(time.time() - start_Time)
"""
