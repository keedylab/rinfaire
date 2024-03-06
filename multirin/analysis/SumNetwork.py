import xarray as xr
import pandas as pd
import networkx as nx
import numpy as np
from pyvis.network import Network
import logging
import pickle
import csv
from multirin.generate.Structure import Structure
from multirin.generate.Subset import generateSubsets

class SumNetwork:

    def __init__ (self, args):
        self.args = args

    def readPickle (self):
        
        # Opens pickle file
        with open(self.args.filename, 'rb') as pickleFile:
            self.multinet = pickle.load(pickleFile)

    def generateSumNetworkAll (self):

        """
        Calculates the sum network for all structures in the MultiNetwork. Main running function for sum network for all structures.
        """

        self.calculateSum(self.multinet.array)
        logging.info(f'Created the sum matrix of all networks')

    def generateSumNetworkSubset (self):

        """
        Calculates the sum network for each subset. Main running function for sum networks of subsets.
        """

        # Generates subsets using the generateSubsets function in Subset.py
        subsetArrays = generateSubsets(self.multinet, self.args.subset)
        print(subsetArrays)

    def calculateSum (self, inputArray):

        """
        This function calculates the sum across all networks when given an input 3D array

        Input: Three dimensional array with axes [network, resi1, resi2]
        Output: Two dimensional array with axes [resi1, resi2] after summing across the network axis
        """

        # Calculates the sum across the network dimension (i.e. for each i,j residue pair)
        sumArray = inputArray.sum(dim="network")

        # Scales the sum array to all be values between 0 and 20
        if self.args.no_scale_sum_network == False:
            sumArray = self.scaleSumNetwork(sumArray)
            logging.info(f'Scaled the Sum Network to values between 0 and 20')

        # Removes weak edges so that the network is easier to visualize
        if self.args.remove_weak_edges != None:
            sumArray = self.removeWeakEdges(sumArray)
            logging.info(f'Removed weak edges from array')

        return sumArray

    def scaleSumNetwork (self, sumArray):

        # Gets maximum value across all dimensions
        maxValue = sumArray.max(dim=['firstResi','secondResi']).item()

        # Then divides each network by max value
        # Scales to a set value (default 0 to 20)
        sumArray = (sumArray / maxValue) * self.args.sum_network_scale

        return sumArray

    def removeWeakEdges (self, sumArray):
        
        # Sorts the array and finds the maximum index to keep (by taking size of array * percent to cutoff)
        # Then gets value at this index
        sortedArray = np.sort(sumArray, axis=None)
        sortedArray = sortedArray[sortedArray != 0]
        maxIndex = round(sortedArray.size * (self.args.remove_weak_edges / 100))
        cutoffValue = sortedArray[maxIndex]

        # Finds entries in the array where they are less than the cutoff threshold
        # Then it replaces them with 0.0
        # If not, then it keeps the original value
        sumArray = xr.where(sumArray < cutoffValue, 0.0, sumArray)

        return sumArray

    def constructGraph (self):

        # Converts XArray into Numpy array
        nparray = self.sumArray.to_numpy()

        # Creates graph from Numpy array
        self.graph = nx.from_numpy_array(nparray)
        self.graph.remove_nodes_from(list(nx.isolates(self.graph)))
        nx.convert_node_labels_to_integers(self.graph)
        
        for i in self.graph.nodes():
            self.graph.nodes[i]['label'] = str(i)

        # Resizing nodes by the degree of the node
        if self.args.no_resize_by_degree == False:
            self.graph = self.resizeByDegree(self.graph)
            
        # Removes subgraphs with < n nodes
        if self.args.remove_subgraphs != 0:
            self.graph = self.removeSubGraphs(self.graph)

        # Shifts sequence to reference sequence
        if self.args.seq_to_ref != None:
            self.graph = self.seqToRef(self.graph)

        # Detects communities within sum graph
        if self.args.detect_communities == True:
            self.graph = self.detectCommunities(self.graph)   

    def visualize (self):    
 
        # Sets PyVis representation
        nts = Network(notebook=True)
        
        # populates the nodes and edges data structures
        nts.from_nx(self.graph)

        # Set deterministic network position using the Kamada-Kawai network layout
        # Solution from: https://stackoverflow.com/questions/74108243/pyvis-is-there-a-way-to-disable-physics-without-losing-graphs-layout
        pos = nx.kamada_kawai_layout(self.graph, scale=2000)

        for node in nts.get_nodes():
            nts.get_node(node)['x']=pos[node][0]
            nts.get_node(node)['y']=-pos[node][1] #the minus is needed here to respect networkx y-axis convention 
            nts.get_node(node)['physics']=False
            nts.get_node(node)['label']=str(node) #set the node label as a string so that it can be displayed
   
        # Outputs the network graph
        nts.toggle_physics(True)
        nts.show_buttons(filter_=['nodes'])
        outputpath = f'{self.args.outputname}.html'
        nts.show(outputpath)

    def resizeByDegree (self, G):

        # Gets the degree of each node and stores this in a dict, then scales it by factor
        # Code from: https://stackoverflow.com/questions/70438752/dynamic-node-sizes-in-pyvis
        nodeDegrees = dict(G.degree(weight='weight')) 
        nodeDegrees.update((x, self.args.resize_by_degree_scale * y) for x, y in nodeDegrees.items())

        # Then sets the node attribute in networkX of size to be the degree
        nx.set_node_attributes(G,nodeDegrees,'size')

        return G
    
    def removeSubGraphs (self, G):

        # Removes components of length less than self.args.remove_subgraphs
        # Code from: https://stackoverflow.com/questions/38308865/how-to-remove-small-components-from-a-graph
        for component in list(nx.connected_components(G)):
            if len(component) < self.args.remove_subgraphs:
                for node in component:
                    G.remove_node(node)

        return G
    
    def seqToRef (self, G):

        # Creates new Structure object from inputted reference structure
        # This is only needed because of the sequenceList attribute in Structure that is neccessary for AllToOne()
        structName = self.args.seq_to_ref
        struct = Structure(structName, self.args)
        
        # Iterates over each node in the graph
        # Then uses the AllToOne function to shift residue numbering back from the alignment positions to the reference sequence numbers
        mappingDict = {}
        for i in G.nodes:

            newLabel = int(self.allToOne(self.multinet.seqaln, struct.name, struct.sequenceList, int(G.nodes[i]['label'])))
            mappingDict[int(G.nodes[i]['label'])] = newLabel
            G.nodes[i]['label'] = newLabel

        # Uses networkX relabel function to relabel nodes using this mapping dictionary
        G = nx.relabel_nodes(G, mappingDict, copy=True)
        return G

    def allToOne (self, seqaln, seqID, sequenceList, mainResidue):

        mainCount = 0
        seqIndex = 0
        
        # Iterates over each residue in the sequence of the structure being queried
        for i in seqaln[seqID]:

            # Increments the main count that counts the total number of positions it has passed on the alignment
            mainCount = mainCount + 1

            # If there is a residue present (not just a - used as a placeholder)
            # Then increment the seq index that counts the total number of actual residue positions that have been passed
            if i != '-':
                seqIndex = seqIndex + 1
                #print(seqIndex)

            # Condition when it reaches the residue in question 
            # (since seq index represents the position in the individual sequence and not the full alignment)
            if mainCount == mainResidue:
            
                # Returns the main count which is the position on the full alignment
                return(sequenceList[seqIndex])
            
    def detectCommunities (self, G):

        # Community detection using the Girvan Newman method, returns a list of sets of communities for every k
        communities = list(nx.community.girvan_newman(G))

        # Modularity -> measures the strength of division of a network into modules
        # Pandas DF that stores modularity for every number of communities (k)
        # From networkX docs: https://networkx.org/documentation/stable/auto_examples/algorithms/plot_girvan_newman.html#sphx-glr-auto-examples-algorithms-plot-girvan-newman-py
        modularity_df = pd.DataFrame(
            [
                [len(communities[k]), nx.community.modularity(G, communities[k])]
                for k in range(len(communities))
            ],
            columns=["k", "modularity"],
        )
        
        # Outputs the modularity df as a csv file if the user specifies
        if self.args.output_modularity == True:
            modularity_df.to_csv(f"{self.args.outputname}_Modularity.csv")

        # Only keeps rows that are smaller than the maximum number of groups (if specified by user)
        if self.args.n_communities != None:
            modularity_df = modularity_df[modularity_df['k'] == self.args.n_communities]

        # Selects the classification with the highest modularity score
        idxModularityMax = modularity_df['modularity'].idxmax()
        modularityMax = modularity_df['modularity'].max()
        selectedCommunities = communities[idxModularityMax]

        # Then labels each node in each community with an associated group
        # Pyvis then colors these groups separately during visualization
        communityCounter = 0

        for community in selectedCommunities:

            nodeList = []
            for node in community:
                G.nodes[node]['group'] = communityCounter
                nodeList.append(int(G.nodes[node]['label']))

            print(f'Community {communityCounter + 1}: {nodeList}')

            communityCounter += 1

        print(f"Modularity of the communities: {modularityMax}")
        
        return G

    def getDegreeInfo (self):

        # Gets degrees of every node and stores it as a dataframe
        nodeDegrees = list(self.graph.degree(weight='weight'))
        nodeDegreesDF = pd.DataFrame.from_records(nodeDegrees, columns = ['Resi,', 'Degree'])
        return nodeDegreesDF
    
    def getLabelInfo (self):

        # Gets label (from clustering) of every node
        nodeLabel = list(self.graph.nodes.data('group'))
        nodeLabelDF = pd.DataFrame.from_records(nodeLabel, columns = ['Resi,', 'Label'])
        return nodeLabelDF

    def getEdgeInfo (self):

        # Gets edge weights
        graphWeights = list(self.graph.edges.data('weight'))
        graphWeightsDF = pd.DataFrame.from_records(graphWeights, columns = ['Resi_1,', 'Resi_2', 'Weight'])
        return graphWeightsDF

    def exportGraphInfo (self):

        import matplotlib.pyplot as plt

        # Creates dataframes for each attribute
        degreeInfo = self.getDegreeInfo()
        edgeInfo = self.getEdgeInfo()
        labelInfo = self.getLabelInfo()

        # Exports these dataframes as CSV files
        filenameDegree = f'{self.args.outputname}_DegreeInfo'
        filenameEdge = f'{self.args.outputname}_EdgeInfo'
        filenameLabel = f'{self.args.outputname}_LabelInfo'

        degreeInfo.to_csv(filenameDegree + '.csv')
        edgeInfo.to_csv(filenameEdge + '.csv')
        labelInfo.to_csv(filenameLabel + '.csv')

        # Plotting histograms for degree and edge weights
        binWidth = 5
        degreeInfo.hist(column='Degree', bins=np.arange(0, degreeInfo['Degree'].max() + binWidth, binWidth))
        plt.savefig(filenameDegree + '.png')

        binWidth = 1
        edgeInfo.hist(column='Weight', bins=np.arange(0, edgeInfo['Weight'].max() + binWidth, binWidth))
        plt.savefig(filenameEdge + '.png')


        # # Writes a csv file of the degree dictionary
        # # Solution from: https://stackoverflow.com/questions/10373247/how-do-i-write-a-python-dictionary-to-a-csv-file
        # filenameCSVDegree = f'{self.args.outputname}_DegreeInfo.csv'

        # with open(filenameCSVDegree,'w') as f:
        #     degreeCSV = csv.writer(f)
        #     degreeCSV.writerows(degreeInfo.items())

        # # Writes a csv file of the edge list
        # # Solution from: https://stackoverflow.com/questions/15578331/save-list-of-ordered-tuples-as-csv
        # filenameCSVEdge = f'{self.args.outputname}_EdgeInfo.csv'    
        
        # with open(filenameCSVEdge,'wb') as out:
        #     edgeCSV = csv.writer(out)
        #     for row in edgeInfo:
        #         edgeCSV.writerow(row)

    def exportPickle (self):

        # Creates new pickle (.pkl) file and then dumps the entire class object into the pickle file
        with open(f'{self.args.outputname}.pkl', 'wb') as pickleFile:
            pickle.dump(self, pickleFile)
