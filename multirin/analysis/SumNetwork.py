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
        self.sumArrays = {}
        self.graphs = {}

    def readPickle (self):
        
        # Opens pickle file
        with open(self.args.filename, 'rb') as pickleFile:
            self.multinet = pickle.load(pickleFile)

    def generateSumNetworkAll (self):

        """
        Calculates the sum network for all structures in the MultiNetwork. Main running function for sum network for all structures.

        Input:
        - self.multinet - object that has MultiNetwork of all structures

        Output:
        - Adds values to self.sumArrays with key 'All'
        """
        
        # Calculates all network sum array and appends it to dictionary of sumArrays with key 'All'
        self.sumArrays['All'] = self.calculateSum(self.multinet.array)
        logging.info(f'Created the sum array of all networks')

    def generateSumNetworkSubset (self, classifier, groupName=None):

        """
        Calculates the sum network for each subset. Main running function for sum networks of subsets.

        Inputs:
        - self.multinet - object that has MultiNetwork of all structures
        - classifier - column of metadata to do subsetting by
        - groupName – name of the specific group in the column that is of interest (optional argument)

        Output:
        - Adds values to self.sumArrays which are all the sum arrays with keys of the group
        """

        # Generates subsets using the generateSubsets function in Subset.py
        subsetMultiNetworks = generateSubsets(self.multinet, classifier, groupName=groupName)

        # Loops over each subset MultiNetwork object in the list
        for subsetMultiNetwork in list(subsetMultiNetworks):

            # Calculates the sum array (takes into account all flags provided)
            self.sumArrays[subsetMultiNetwork] = self.calculateSum(subsetMultiNetworks[subsetMultiNetwork].array)
            
            # Tests whether the output from calculateSum is None (meaning the sumArray was empty)
            if self.sumArrays[subsetMultiNetwork] is None:

                # If so, then delete MultiNetwork and sumArray related to that group
                print(f'Sum Array for {subsetMultiNetwork} is empty, removing from dictionary')
                del subsetMultiNetworks[subsetMultiNetwork]
                del self.sumArrays[subsetMultiNetwork]

            else:
                logging.info(f'Created the subset sum array of {subsetMultiNetwork}')

    def calculateSum (self, inputArray):

        """
        This function calculates the sum across all networks when given an input 3D array

        Input: Three dimensional array with axes [network, resi1, resi2]
        Output: Two dimensional array with axes [resi1, resi2] after summing across the network axis
        """

        # Calculates the sum across the network dimension (i.e. for each i,j residue pair)
        sumArray = inputArray.sum(dim="network")

        # Tests if the sum of the entire sumArray = 0 (array is empty with no values)
        # If so, function ends and returns None
        if sumArray.sum() == 0:
            return None

        # Scales the sum array to all be values between 0 and 20
        # if self.args.no_scale_sum_network == False:
        if self.args.scale_sum_network == 'max':
            sumArray = self.scaleSumNetworkMax(sumArray)
            logging.info(f'Scaled the Sum Network to values between 0 and 20')
        elif self.args.scale_sum_network == 'struct':
            sumArray = self.scaleSumNetworkStruct(sumArray, inputArray)
            logging.info(f'Scaled the Sum Network by total number of structures')
        elif self.args.scale_sum_network == 'none':
            pass
        else:
            raise NameError("Scaling type not found")

        # Removes weak edges so that the network is easier to visualize
        if self.args.remove_weak_edges != None:
            sumArray = self.removeWeakEdges(sumArray)
            logging.info(f'Removed weak edges from array')

        # Confirming the sum array is symmetric with diagonal entries = 0
        sumNPTest = sumArray.to_numpy()
        if ((np.all(np.diag(sumNPTest)==0)) or (np.allclose(sumNPTest.T, sumNPTest))) == False:
            print('Array is not symmetric with diagonal entries of zero')
            return None

        return sumArray

    def scaleSumNetworkMax (self, sumArray):

        """
        Scales the Sum Network by the maximum value so max value is the scaling factor
        """

        # Gets maximum value across all dimensions
        maxValue = sumArray.max(dim=['firstResi','secondResi']).item()

        # Then divides each network by max value
        # Scales to a set value (default 0 to 20)
        sumArray = (sumArray / maxValue) * self.args.sum_network_scaling_factor

        return sumArray
    
    def scaleSumNetworkStruct (self, sumArray, originalArray):

        """
        Scales the Sum Network by the number of total structures that were summed
        """

        # Gets total number of structures/networks
        numberOfStructures = originalArray.sizes['network']
        
        # Then divides each network by total number of structures, also scales by scaling factor
        sumArray = (sumArray / numberOfStructures) * self.args.sum_network_scaling_factor

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
    
    def constructMaxSpanningTrees (self):

        """
        Function that creates the maximum spanning tree as a way to visualize the sum network. Done instead of the traditional method of visualizing the sum network that shows all edges.

        Inputs:
        - self.sumArrays: Dictionary of sum arrays for each subset

        Outputs:
        - self.graphs: Dictionary with entries being the graphs for each subset
        """

        # Loops through each sum array across all subsets (if no subsets it just iterates once on the full sum array)
        for sumArray in self.sumArrays:

            # Converts XArray into Numpy array
            nparray = self.sumArrays[sumArray].to_numpy()

            # Creates graph from Numpy array
            newGraph = nx.from_numpy_array(nparray)
            newGraph.remove_nodes_from(list(nx.isolates(newGraph)))

            for i in newGraph.nodes():
                newGraph.nodes[i]['label'] = str(i)

            # Finds maximum spanning tree
            MSTGraph = nx.maximum_spanning_tree(newGraph)

            # Resizing nodes by the degree of the node
            if self.args.no_resize_by_degree == False:
                MSTGraph = self.resizeByDegree(MSTGraph)

            # Shifts sequence to reference sequence
            if self.args.seq_to_ref != None:
                MSTGraph = self.seqToRef(MSTGraph)

            # Detects communities within sum graph
            if self.args.detect_communities == True:
                MSTGraph = self.detectCommunities(MSTGraph) 

            # Appends this new graph to the dictionary of graphs for each sum array, with the key being the original name in the sumArrays dictionary
            self.graphs[sumArray] = MSTGraph

    def constructGraphs (self):

        """
        Function that generates networkx graphs for each sum array.

        Input: self.sumArrays dictionary with all the sumArrays
        Output: New dictionary entries in self.graphs which are all the sum graphs
        """

        # Loops over each sum array (in the event there's > 1 sum array in the case of subsets)
        for sumArray in self.sumArrays:

            # Converts XArray into Numpy array
            nparray = self.sumArrays[sumArray].to_numpy()

            # Creates graph from Numpy array
            newGraph = nx.from_numpy_array(nparray)
            newGraph.remove_nodes_from(list(nx.isolates(newGraph)))
            
            for i in newGraph.nodes():
                newGraph.nodes[i]['label'] = str(i)

            # Resizing nodes by the degree of the node
            if self.args.no_resize_by_degree == False:
                newGraph = self.resizeByDegree(newGraph)
                
            # Removes subgraphs with < n nodes
            if self.args.remove_subgraphs != 0:
                newGraph = self.removeSubGraphs(newGraph)

            # Shifts sequence to reference sequence
            if self.args.seq_to_ref != None:
                newGraph = self.seqToRef(newGraph)

            # Detects communities within sum graph
            if self.args.detect_communities == True:
                newGraph = self.detectCommunities(newGraph)   

            # Appends this new graph to the dictionary of graphs for each sum array, with the key being the original name in the sumArrays dictionary
            self.graphs[sumArray] = newGraph

    def visualize (self):

        """
        Function that generates the PyVis visualization of each networkx graph.
        
        Input: self.graphs dictionary of networkx graphs
        Output: .html file visualizations of each network
        """

        # Iterates over each networkx sum graph
        for graph in self.graphs:    
 
            # Sets PyVis representation
            nts = Network(notebook=True)
            
            # populates the nodes and edges data structures
            nts.from_nx(self.graphs[graph])

            # Set deterministic network position using the Kamada-Kawai network layout
            # Solution from: https://stackoverflow.com/questions/74108243/pyvis-is-there-a-way-to-disable-physics-without-losing-graphs-layout
            pos = nx.kamada_kawai_layout(self.graphs[graph], scale=2000)

            for node in nts.get_nodes():
                nts.get_node(node)['x']=pos[node][0]
                nts.get_node(node)['y']=-pos[node][1] #the minus is needed here to respect networkx y-axis convention 
                # nts.get_node(node)['physics']=False
                nts.get_node(node)['label']=str(node) #set the node label as a string so that it can be displayed
    
            # Outputs the network graph
            nts.toggle_physics(True)
            nts.show_buttons(filter_=['nodes'])
            outputpath = f'{self.args.outputname}_{graph}.html'
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

            newLabel = self.allToOne(self.multinet.seqaln, struct.name, struct.sequenceList, int(G.nodes[i]['label']))
            
            mappingDict[int(G.nodes[i]['label'])] = newLabel
            G.nodes[i]['label'] = newLabel

        # Uses networkX relabel function to relabel nodes using this mapping dictionary
        G = nx.relabel_nodes(G, mappingDict, copy=True)

        # Removes nodes that didn't map to residue in structure
        if self.args.keep_nan == False:
            remove = [node for node in G.nodes if str(node)[0:3] == 'NaN']
            G.remove_nodes_from(remove)

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

                if i == '-':
                    return(f'NaN_{str(sequenceList[seqIndex])}_{str(mainCount)}')

                else:
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
        communityList = []
        nodeList = []

        # Iterates over each community
        for community in selectedCommunities:

            # Iterates over each residue in community
            communityNodeList = []
            for node in community:

                # Sets the nodes' group attribute to the be the community number
                G.nodes[node]['group'] = communityCounter

                # Adds to lists to output
                communityNodeList.append(int(G.nodes[node]['label']))
                nodeList.append(node)
                communityList.append(communityCounter + 1)

            # Prints community with list of residues in it plus adds it to dictionary
            print(f'Community {communityCounter + 1}: {communityNodeList}')

            # updates community number
            communityCounter += 1

        print(f"Modularity of the communities: {modularityMax}")

        # Output communities as csv file
        communityDF = pd.DataFrame.from_dict({'Residue': nodeList, 'Community':communityList}).set_index('Residue')
        communityDF.to_csv(f'{self.args.outputname}Communities.csv')
        
        return G

    def getDegreeInfo (self, graph):

        # Gets degrees of every node and stores it as a dataframe
        nodeDegrees = list(graph.degree(weight='weight'))
        nodeDegreesDF = pd.DataFrame.from_records(nodeDegrees, columns = ['Resi,', 'Degree'])
        return nodeDegreesDF
    
    def getLabelInfo (self, graph):

        # Gets label (from clustering) of every node
        nodeLabel = list(graph.nodes.data('group'))
        nodeLabelDF = pd.DataFrame.from_records(nodeLabel, columns = ['Resi,', 'Label'])
        return nodeLabelDF

    def getEdgeInfo (self, graph):

        # Gets edge weights
        graphWeights = list(graph.edges.data('weight'))
        graphWeightsDF = pd.DataFrame.from_records(graphWeights, columns = ['Resi_1,', 'Resi_2', 'Weight'])
        return graphWeightsDF

    def exportGraphInfo (self):

        """
        This function exports information about the sum network graphs (such as node degree, edges, and labels)
        """

        import matplotlib.pyplot as plt

        # Loops over all graphs
        for graph in self.graphs:

            # Creates dataframes for each attribute
            degreeInfo = self.getDegreeInfo(self.graphs[graph])
            edgeInfo = self.getEdgeInfo(self.graphs[graph])
            labelInfo = self.getLabelInfo(self.graphs[graph])

            # Exports these dataframes as CSV files
            filenameDegree = f'{self.args.outputname}_{graph}_DegreeInfo'
            filenameEdge = f'{self.args.outputname}_{graph}_EdgeInfo'
            filenameLabel = f'{self.args.outputname}_{graph}_LabelInfo'

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

    def exportPickle (self):

        # Creates new pickle (.pkl) file and then dumps the entire class object into the pickle file
        with open(f'{self.args.outputname}.pkl', 'wb') as pickleFile:
            pickle.dump(self, pickleFile)
