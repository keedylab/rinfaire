import networkx as nx
import numpy as np
import pandas as pd
import gemmi
from pyvis.network import Network
import logging
import pickle
import csv
from multirin.generate.Structure import Structure
from multirin.generate.IndividualNetwork import IndividualNetwork
from argparse import Namespace

class ResiduesOfInterest:

    def __init__ (self, args):
        self.args = args

    def readPickle (self):
        
        # Opens pickle file
        with open(self.args.filename, 'rb') as pickleFile:
            self.sumNetwork = pickle.load(pickleFile)

    def readInputSetFile (self):

        # Opens .csv file as pandas dataframe
        dfInputSet = pd.read_csv(self.args.input_set)

        # Creates dict of lists to store every set of residues
        self.inputSetDict = {}

        for col in dfInputSet.columns:
            # Gets column, drops N/A values, and converts values to ints
            setColumn = dfInputSet[col]
            setColumn = setColumn.dropna()
            setColumn = setColumn.astype(int)

            # Converts to list and adds to the dictionary
            self.inputSetDict[col] = setColumn.to_list()

    def findOverlapInputSet (self):

        # Creates a dictionary to store intersecting residues
        self.overlapDict = {}

        # Condition if to include adjacent residues to the network or not
        if self.args.include_adjacent_residues != None:

            # Sets input structure file as Structure object
            inputStruct = Structure(self.args.include_adjacent_residues, None)

            # Creates IndividualNetwork object using the sumNetwork as an input network and inputStruct from user as reference structure
            # Then uses the addAdjacentResidues() method to find adjacent residues to the sumNetwork
            args = Namespace(no_norm_resi=False)
            self.adjResisNetwork = IndividualNetwork(inputStruct, args, network=self.sumNetwork.graph)
            self.adjResisNetwork.addAdjacentResidues()

            networkList = list(self.adjResisNetwork.network.nodes)

        else:
            networkList = list(self.sumNetwork.graph.nodes)

        # Iterates over each sector of residues in the input set dictionary
        for col in self.inputSetDict:

            # Intersection between two lists
            intersectionList = [value for value in self.inputSetDict[col] if value in networkList]

            # Finds the percent overlap between the intersection and the total length of the input set
            overlapPercent = (len(intersectionList) / len(self.inputSetDict[col])) * 100

            # Prints stats out
            addString = ""
            if self.args.include_adjacent_residues != None:
                addString = "(and adjacent to)"

            print(f'{col}: \n   {overlapPercent}% of residues are found in {addString} network \n   Common residues are: {intersectionList} \n')

            # Appends list to overlap dictionary
            self.overlapDict[col] = intersectionList

    def labelGraphOverlap (self):

        # Creates a copy of the graph to plot the overlapping residues
        self.overlapGraph = self.sumNetwork.graph

        # Then labels each node in each community with an associated group
        # Pyvis then colors these groups separately during visualization
        communityColors = ['#0077BB', '#EE7733', '#33BBEE', '#EE3377', '#CC3311', '#009988']
        communityCounter = 0

        # Sets all node colors to gray by default
        for node in self.overlapGraph.nodes:
            self.overlapGraph.nodes[node]['color'] = '#BBBBBB'

        # Then recolors them by community
        for community in self.overlapDict:

            for node in self.overlapDict[community]:
                self.overlapGraph.nodes[node]['color'] = communityColors[communityCounter]

            communityCounter += 1

    def visualize (self, graph, filename):    
 
        # Sets PyVis representation
        nts = Network(notebook=True, width="100%", height="50vw")
        nts.set_options("""
        var options = {
        "nodes": {
            "font": {
            "size": 25,
            "face": "arial",
            "physics": false
            }
        }
        }
        """)
        
        # populates the nodes and edges data structures
        nts.from_nx(graph)

        # Set deterministic network position using the Kamada-Kawai network layout
        # Solution from: https://stackoverflow.com/questions/74108243/pyvis-is-there-a-way-to-disable-physics-without-losing-graphs-layout
        pos = nx.kamada_kawai_layout(graph, scale=2000)

        for node in nts.get_nodes():
            nts.get_node(node)['x']=pos[node][0]
            nts.get_node(node)['y']=-pos[node][1] #the minus is needed here to respect networkx y-axis convention 
            nts.get_node(node)['physics']=False
            nts.get_node(node)['label']=str(node) #set the node label as a string so that it can be displayed
   
        # Outputs the network graph
        outputpath = f'{self.args.outputname}_{filename}.html'
        nts.show(outputpath)

    def exportPickle (self):

        # Creates new pickle (.pkl) file and then dumps the entire class object into the pickle file
        with open(f'{self.args.outputname}.pkl', 'wb') as pickleFile:
            pickle.dump(self, pickleFile)
