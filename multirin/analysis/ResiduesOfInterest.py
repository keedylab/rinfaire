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
import random

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

        networkList = list(self.sumNetwork.nodes)

        # Iterates over each sector of residues in the input set dictionary
        for col in self.inputSetDict:

            # Intersection between two lists
            inputSet = set(self.inputSetDict[col])
            networkSet = set(networkList)
            intersectionSet = inputSet.intersection(networkSet)
            intersectionList = sorted(list(intersectionSet))

            # Find values not in each respective set
            networkNotInputList = sorted(list(networkSet.difference(intersectionSet)))
            inputNotNetworkList = sorted(list(inputSet.difference(intersectionSet)))
            
            #intersectionList = [value for value in self.inputSetDict[col] if value in networkList]

            # Finds the percent overlap between the intersection and the total length of the input set
            overlapPercent = (len(intersectionList) / len(self.inputSetDict[col])) * 100

            # Prints stats out
            addString = ""
            if self.args.include_adjacent_residues != None:
                addString = "(and adjacent to)"

            print(f'''
{col}: \n   
    {overlapPercent}% of residues are found in {addString} network \n   
    Has {len(intersectionList)} residues in common: {intersectionList} \n
    Has {len(networkNotInputList)} residues in network but not input set: {networkNotInputList} \n
    Has {len(inputNotNetworkList)} residues in input set but not network: {inputNotNetworkList} \n
            ''')

            # Appends list to overlap dictionary
            self.overlapDict[col] = intersectionList

    def findSignificance (self):

        import scipy
        import statistics

        ### First gets a network of all possible contacts between residues in input structure
        # Sets input structure file as Structure object
        inputStruct = Structure(self.args.find_significance, None)

        # Creates IndividualNetwork object using the sumNetwork as an input network and inputStruct from user as reference structure
        # Then uses the addallResidues() method to find all possible residue - residue contacts
        args = Namespace(no_norm_resi=False)
        self.allResisNetwork = IndividualNetwork(inputStruct, args)
        self.allResisNetwork.addAllResidues()

        allNetworkList = list(self.allResisNetwork.network.nodes)

        ### Then finds the fraction of residues to each input set residue that are close to network residues
        closeInputResiCountList = self.findFractionCloseToNetwork(self.inputSetDict[self.args.col])

        # Finds the set of residues that are not in the input set
        #nonInputResiList = [i for i in allNetworkList if i not in self.inputSetDict[self.args.col]]

        # ### Then finds the fraction of residues to each non input set residue that are close to network residues
        # closeNonInputResiCountList = self.findFractionCloseToNetwork(nonInputResiList)

        sig_test_values = []

        for i in range(0, int(self.args.n_iter_sig_test)):

            # Finds random sample of residues
            #print(len(self.inputSetDict[self.args.col]))
            randomResiList = random.sample(allNetworkList, len(self.inputSetDict[self.args.col])) # allNetworkList
            #print(randomResiList)

            ### Then finds the fraction of residues to each non input set residue that are close to network residues
            randomResiCountList = self.findFractionCloseToNetwork(randomResiList)

            # Finds the statistical significance between all input set vs non input set residues' fraction of residues close to network
            # Does this by a KS test
            ksStat, pValue = scipy.stats.ks_2samp(closeInputResiCountList, randomResiCountList) #closeNonInputResiCountList

            sig_test_values.append({'input_set_mean': statistics.mean(closeInputResiCountList), 
                                    'non-input_set_mean': statistics.mean(randomResiCountList), 
                                    'ks_stat': ksStat,
                                    'pvalue': pValue})

            # print(f'Input set mean: {statistics.mean(closeInputResiCountList)}')
            # print(f'Non-input set mean: {statistics.mean(randomResiCountList)}')
            # print(f'The p value is: {pValue}')

        sig_test_df = pd.DataFrame(sig_test_values)
        sig_test_df.to_csv(f'{self.args.outputname}_{self.args.col}_sig_test.csv', index=False)

        # Plots histograms if either flag is provided
        if (self.args.cumulative_histogram or self.args.histogram) == True:
            self.plotHistogram(closeInputResiCountList, randomResiCountList)
            self.plotBoxPlot(closeInputResiCountList, randomResiCountList)
    
    def plotHistogram (self, closeInputResiCountList, closeNonInputResiCountList):

        """
        This function plots the histograms between the input set given and non input set
        It is with respect to the fraction/number of residues that are adjacent to the network in the input vs non-input set

        Inputs:
        - closeInputResiCountList: The list of fractions/numbers of the number of adjacent residues that are in the network for each residue in the input set
        - closeNonInputResiCountList: Same as closeInputResiCountList but for the non-input set

        Outputs:
        - Matplotlib plots of both histograms as .png files
        """

        import matplotlib.pyplot as plt

        # # Now plots the two distributions as histograms
        # fig, ax1 = plt.subplots()

        # Sets the bin range for the histogram
        if self.args.no_normalize_by_total == False:
            binRange = np.arange(0.0, 1.0, 0.05)
        
        else:
            maxData = max([max(closeInputResiCountList),max(closeNonInputResiCountList)])
            binWidth = 1
            binRange = np.arange(0, maxData+binWidth, binWidth)
    
        # Creates two histogram plots
        color = 'blue'
        plt.hist(closeInputResiCountList, label=f'{self.args.col} Residues', color=color, bins=binRange, density=True, cumulative=self.args.cumulative_histogram, histtype='step', linewidth=1.5)
        
        color = 'red'
        plt.hist(closeNonInputResiCountList, label=f'Non {self.args.col} Residues', color=color, bins=binRange, density=True, cumulative=self.args.cumulative_histogram, histtype='step', linewidth=1.5)

        # Adds the axis labels and figure legend
        if self.args.no_normalize_by_total == False:
            plt.xlabel("Fraction of adjacent residues (within 4Å) in network")
    
        else:
            plt.xlabel("Number of adjacent residues (within 4Å) in network")

        if self.args.cumulative_histogram == True:
            plt.ylabel("Cumulative Fraction")
            plt.legend(loc="center right")

        else:
            plt.ylabel("Total Fraction")
            plt.legend(loc="upper right")

        # Saving figure
        plt.title(self.args.col)
        plt.savefig(f'{self.args.outputname}_{self.args.col}_Histogram.png')

    def plotBoxPlot (self, closeInputResiCountList, closeNonInputResiCountList):

        import matplotlib.pyplot as plt
        import numpy as np

        # Create a figure and axes
        fig, ax = plt.subplots(figsize=(4, 6))

        # Create the boxplots with customized colors
        boxprops = dict(color='black', linewidth=0.5)  # Gray outline for boxes
        whiskerprops = dict(color='black', linewidth=0.5)  # Whisker thickness
        capprops = dict(color='black', linewidth=0.5)  # Cap thickness
        medianprops = dict(color='black', linewidth=0.5)  # Gray median line

        # Create the boxplots
        ax.boxplot([closeInputResiCountList, closeNonInputResiCountList], labels=[f'{self.args.col}', f'Random Set'], widths=0.5, boxprops=boxprops, medianprops=medianprops, whiskerprops=whiskerprops, capprops=capprops)

        # Add title and labels
        plt.xlabel('Residues', fontsize=15)
        plt.ylabel('Network Residues within 4Å', fontsize=15)

        # Get min/max for y-axis ticks
        all_values = closeInputResiCountList + closeNonInputResiCountList
        min_value, max_value = min(all_values), max(all_values)

        # Set y-axis ticks at step of 2
        plt.yticks(np.arange(min_value, max_value + 2, 2), fontsize=15)
        plt.xticks(fontsize=15)
        combinedLists = [closeInputResiCountList, closeNonInputResiCountList]
        colors = ['darkorange', 'darkgray']


        for i, (y, color) in enumerate(zip(combinedLists, colors)):
            y = combinedLists[i]
            x = np.random.normal(1+i, 0.04, size=len(y))
            plt.plot(x, y, '.', color=color, alpha=0.5)
            plt.tight_layout()

        # Show the plot
        plt.savefig(f'{self.args.outputname}_{self.args.col}_BoxPlot.png', bbox_inches='tight', dpi=300)

    def findFractionCloseToNetwork (self, resiList):
        
        closeResiCountList = []
        
        # Loops over every residue of interest in the input set
        for resi in resiList:
            
            # Checks if the graph has the node or not
            # If the node is not present then it is not close to any of the network residues
            if self.allResisNetwork.network.has_node(resi) == False:
                closeResiCountList.append(0)

            # Otherwise, it has at least one connection, and the number of connections is the degree
            else:

                # Finds the degree of the adjacent residue network (to find number of network adjacent connections)
                closeNetworkResiCount = self.conditionalDegree(resi)
                
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

                # print(resi, closeNetworkResiCount)

        return(closeResiCountList)
    
    def conditionalDegree (self, inputNode):

        """
        Function that gets the degree of a node where its edges must pass a certain condtion
        Ex. get only degree of node to only network connections (network residues)

        Based off of solution from: https://stackoverflow.com/questions/30077957/calculate-the-degree-of-nodes-only-including-edges-with-a-specific-attribute-in
        """

        degree = 0
        for u,v,d in self.allResisNetwork.network.edges(inputNode, data=True):

            if ((u == inputNode) & (v in self.sumNetwork.nodes)):
                degree += 1
              
        return degree

    def labelGraphOverlap (self):

        # Creates a copy of the graph to plot the overlapping residues
        self.overlapGraph = self.sumNetwork

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
