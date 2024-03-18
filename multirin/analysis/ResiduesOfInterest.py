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

    def findSignificance (self):

        import scipy
        import statistics
        import matplotlib.pyplot as plt

        ### First gets a network of all possible contacts between residues in input structure
        # Sets input structure file as Structure object
        inputStruct = Structure(self.args.find_significance, None)

        # Creates IndividualNetwork object using the sumNetwork as an input network and inputStruct from user as reference structure
        # Then uses the addallResidues() method to find all possible residue - residue contacts
        args = Namespace(no_norm_resi=False)
        self.allResisNetwork = IndividualNetwork(inputStruct, args)
        self.allResisNetwork.addAllResidues()

        allNetworkList = list(self.allResisNetwork.network.nodes)

        ### Then finds adjacent network
        # Creates IndividualNetwork object using the sumNetwork as an input network and inputStruct from user as reference structure
        # Then uses the addAdjacentResidues() method to find adjacent residues to the sumNetwork
        self.adjResisNetwork = IndividualNetwork(inputStruct, args, network=self.sumNetwork.graph)
        self.adjResisNetwork.addAdjacentResidues()

        ### Then finds the fraction of residues to each input set residue that are close to network residues
        closeInputResiCountList = self.findFractionCloseToNetwork(self.inputSetDict[self.args.col])

        # Finds the set of residues that are not in the input set
        nonInputResiList = [i for i in allNetworkList if i not in self.inputSetDict[self.args.col]]

        ### Then finds the fraction of residues to each non input set residue that are close to network residues
        closeNonInputResiCountList = self.findFractionCloseToNetwork(nonInputResiList)

        # Finds the statistical significance between all input set vs non input set residues' fraction of residues close to network
        # Does this by a student's t-test
        tStat, pValue = scipy.stats.ttest_ind(closeInputResiCountList, closeNonInputResiCountList)

        print(f'Input set mean: {statistics.mean(closeInputResiCountList)}')
        print(f'Non-input set mean: {statistics.mean(closeNonInputResiCountList)}')
        print(f'The p value is: {pValue}')

        # Now plots the two distributions as histograms
        fig, ax1 = plt.subplots()

        # Sets the bin range for the histogram
        if self.args.no_normalize_by_total == False:
            binRange = np.arange(0.0, 1.0, 0.05)
        
        else:
            maxData = max([max(closeInputResiCountList),max(closeNonInputResiCountList)])
            binWidth = 1
            binRange = np.arange(0, maxData+binWidth, binWidth)
    
        color = 'blue'
        ax1.set_ylabel('Input Set Count', color=color)
        ax1.hist(closeInputResiCountList, label='Input Set', alpha=0.5, color=color, bins=binRange)

        ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

        color = 'red'
        ax2.set_ylabel('Non-Input Set Count', color=color)
        ax2.hist(closeNonInputResiCountList, label='Non-Input Set', alpha=0.5, color=color, bins=binRange)

        plt.title(self.args.col)
        plt.savefig(f'{self.args.outputname}_{self.args.col}_Histogram.png')

    def findFractionCloseToNetwork (self, resiList):
        
        closeResiCountList = []
        
        # Loops over every residue of interest in the input set
        for resi in resiList:
            
            # Checks if the graph has the node or not
            # If the node is not present then it is not close to any of the network residues
            if self.adjResisNetwork.network.has_node(resi) == False:
                closeResiCountList.append(0)

            # Otherwise, it has at least one connection, and the number of connections is the degree
            else:

                # Finds the degree of the adjacent residue network (to find number of network adjacent connections)
                closeNetworkResiCount = self.conditionalDegree(self.adjResisNetwork.network, resi, 'edgeClass', None)
                
                # By default, this is true, telling the program to normalize the count of adjacent network residues by the total number of residues
                if self.args.no_normalize_by_total == False:

                    # Then finds the degree of the total residue network (to find number of total connections)
                    closeTotalResiCount = self.allResisNetwork.network.degree[resi]

                    # Divides network adjacent count by total count, appends this to list
                    closeNormCount = closeNetworkResiCount / closeTotalResiCount
                    closeResiCountList.append(closeNormCount)
                
                # Otherwise just add in the number of adjacent network residues
                else:
                    closeResiCountList.append(closeNetworkResiCount)

        return(closeResiCountList)
    
    def conditionalDegree (self, G, inputNode, condition, conditionValue):

        """
        Function that gets the degree of a node where its edges must pass a certain condtion
        Ex. get only degree of node to only network connections (network residues)

        Based off of solution from: https://stackoverflow.com/questions/30077957/calculate-the-degree-of-nodes-only-including-edges-with-a-specific-attribute-in
        """

        degree = 0
        for u,v,d in G.edges(inputNode, data=True):
            if d[condition] == conditionValue:
                degree += 1
                
        return degree

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
