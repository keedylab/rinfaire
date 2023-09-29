import xarray as xr
import networkx as nx
from pyvis.network import Network
import logging
import pickle

class SumNetwork:

    def __init__ (self, args):
        self.args = args
        self.readPickle()

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
        
        # Get maximum value across the entire array
        # Then multiply that by the percent cutoff threshold specified by user
        percentCutoffValue = self.sumArray.max() * (self.args.remove_weak_edges / 100)

        # Finds entries in the array where they are less than the percent cutoff threshold
        # Then it replaces them with 0.0
        # If not, then it keeps the original value
        self.sumArray = xr.where(self.sumArray < percentCutoffValue, 0.0, self.sumArray)

    def visualize (self):

        # Converts XArray into Numpy array
        nparray = self.sumArray.to_numpy()

        # Creates graph from Numpy array
        G = nx.from_numpy_array(nparray)
        G.remove_nodes_from(list(nx.isolates(G)))
        nx.convert_node_labels_to_integers(G)
        
        for i in G.nodes():
            G.nodes[i]['label'] = str(i)
        
        # widths = nx.get_edge_attributes(G, 'weight')

        nts = Network(notebook=True)
        
        # populates the nodes and edges data structures
        nts.from_nx(G)
        outputpath = f'{self.args.outputname}.html'
        nts.show(outputpath)

    def exportPickle (self):

        # Creates new pickle (.pkl) file and then dumps the entire class object into the pickle file
        with open(f'{self.args.outputname}.pkl', 'wb') as pickleFile:
            pickle.dump(self, pickleFile)
