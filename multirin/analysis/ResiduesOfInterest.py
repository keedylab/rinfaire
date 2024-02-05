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

    def exportPickle (self):

        # Creates new pickle (.pkl) file and then dumps the entire class object into the pickle file
        with open(f'{self.args.outputname}.pkl', 'wb') as pickleFile:
            pickle.dump(self, pickleFile)
