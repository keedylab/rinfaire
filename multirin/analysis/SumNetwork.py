import xarray as xr
import networkx as nx
import numpy as np
from pyvis.network import Network
import logging
import pickle

class SumNetwork:

    def __init__ (self, args):
        self.args = args

    def readPickle (self):
        
        # Opens pickle file
        with open(self.args.filename, 'rb') as pickleFile:
            self.multinet = pickle.load(pickleFile)

    def calculateSum (self):

        # Calculates the sum across the network dimension (i.e. for each i,j residue pair)
        self.sumArray = self.multinet.array.sum(dim="network")

        # Scales the sum array to all be values between 0 and 20
        if self.args.no_scale_sum_network == False:
            self.scaleSumNetwork()
            logging.info(f'Scaled the Sum Network to values between 0 and 20')

        # Removes weak edges so that the network is easier to visualize
        if self.args.remove_weak_edges != None:
            self.removeWeakEdges()
            logging.info(f'Removed weak edges from array')

        logging.info(f'Created the sum matrix of all networks')

    def scaleSumNetwork (self):

        # Gets maximum value across all dimensions
        maxValue = self.sumArray.max(dim=['firstResi','secondResi']).item()

        # Then divides each network by max value
        # Scales to a set value (default 0 to 20)
        self.sumArray = (self.sumArray / maxValue) * self.args.sum_network_scale

    def removeWeakEdges (self):
        
        # Sorts the array and finds the maximum index to keep (by taking size of array * percent to cutoff)
        # Then gets value at this index
        sortedArray = np.sort(self.sumArray, axis=None)
        sortedArray = sortedArray[sortedArray != 0]
        maxIndex = round(sortedArray.size * (self.args.remove_weak_edges / 100))
        cutoffValue = sortedArray[maxIndex]

        # Finds entries in the array where they are less than the cutoff threshold
        # Then it replaces them with 0.0
        # If not, then it keeps the original value
        self.sumArray = xr.where(self.sumArray < cutoffValue, 0.0, self.sumArray)

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

    def visualize (self):    

        nts = Network(notebook=True)
        
        # populates the nodes and edges data structures
        nts.from_nx(self.graph)
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

    def exportPickle (self):

        # Creates new pickle (.pkl) file and then dumps the entire class object into the pickle file
        with open(f'{self.args.outputname}.pkl', 'wb') as pickleFile:
            pickle.dump(self, pickleFile)
