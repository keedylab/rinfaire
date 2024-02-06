import networkx as nx
import numpy as np
import pandas as pd
from pyvis.network import Network
import logging
import pickle
import csv

class ResiduesOfInterest:

    def __init__ (self, args):
        self.args = args

    def readPickle (self):
        
        # Opens pickle file
        with open(self.args.filename, 'rb') as pickleFile:
            self.sumNetwork = pickle.load(pickleFile)

    def findOverlapInputSet (self):

        # Opens .csv file as pandas dataframe
        dfInputSet = pd.read_csv(self.args.input_set)

        # Creates a dictionary to store intersecting residues
        self.overlapDict = {}

        # List of nodes from the network
        networkList = list(self.sumNetwork.graph.nodes)

        for col in dfInputSet.columns:
            
            # Gets column, drops N/A values, and converts values to ints
            setColumn = dfInputSet[col]
            setColumn = setColumn.dropna()
            setColumn = setColumn.astype(int)

            # Converts to list
            inputSetList = setColumn.to_list()

            # Intersection between two lists
            intersectionList = [value for value in inputSetList if value in networkList]

            # Finds the percent overlap between the intersection and the total length of the input set
            overlapPercent = (len(intersectionList) / len(inputSetList)) * 100
            print(f'{col}: \n   {overlapPercent}% of residues are found in network \n   Common residues are: {intersectionList} \n')

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

    def visualize (self, graph):    
 
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
   
        for edge in nts.get_edges():
            edge["color"] = '#BBBBBB'

        # Outputs the network graph
        outputpath = f'{self.args.outputname}.html'
        nts.show(outputpath)

    def exportPickle (self):

        # Creates new pickle (.pkl) file and then dumps the entire class object into the pickle file
        with open(f'{self.args.outputname}.pkl', 'wb') as pickleFile:
            pickle.dump(self, pickleFile)
