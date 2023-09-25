import xarray as xr
import networkx as nx
from pyvis.network import Network
import logging

def sum (self):

    # Calculates the sum across the network dimension (i.e. for each i,j residue pair)
    sumArray = self.array.sum(dim="network")

    # Scales the sum array to all be values between 0 and 20
    if self.args.no_scale_sum_network == False:
        sumArray = self.scaleSumNetwork(sumArray)
        logging.info(f'Scaled the Sum Network to values between 0 and 20')

    logging.info(f'Created the sum matrix of all networks')
    return sumArray

def scaleSumNetwork (self, inputArray):

    # Gets maximum value across all dimensions
    maxValue = inputArray.max(dim=['firstResi','secondResi']).item()

    # Then divides each network by max value
    # Scales to a set value (default 0 to 20)
    scaledInputArray = (inputArray / maxValue) * self.args.sum_network_scale

    return scaledInputArray

def removeWeakEdges (self, inputArray):
    
    # Get maximum value across the entire array
    # Then multiply that by the percent cutoff threshold specified by user
    percentCutoffValue = inputArray.max() * (self.args.remove_weak_edges / 100)

    # Finds entries in the array where they are less than the percent cutoff threshold
    # Then it replaces them with 0.0
    # If not, then it keeps the original value
    return xr.where(inputArray < percentCutoffValue, 0.0, inputArray)

def visualize (self, inputarray, outputname):

    # Removes weak edges so that the network is easier to visualize
    if self.args.remove_weak_edges != None:
        inputarray = self.removeWeakEdges(inputarray)
        logging.info(f'Removed weak edges from array for visualization')

    # Converts XArray into Numpy array
    nparray = inputarray.to_numpy()

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
    outputpath = f'{self.args.output}{outputname}.html'
    nts.show(outputpath)