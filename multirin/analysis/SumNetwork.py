import xarray as xr
import networkx as nx
from pyvis.network import Network
import logging
import pickle
from multirin.generate.Structure import Structure

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

        # Resizing nodes by the degree of the node
        if self.args.no_resize_by_degree == False:
            G = self.resizeByDegree(G)    

        # Removes subgraphs with < n nodes
        if self.args.remove_subgraphs != 0:
            G = self.removeSubGraphs(G)

        if self.args.seq_to_ref != None:
            G = self.seqToRef(G)

        # widths = nx.get_edge_attributes(G, 'weight')

        nts = Network(notebook=True)
        
        # populates the nodes and edges data structures
        nts.from_nx(G)
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
        for i in G.nodes:
            G.nodes[i]['label'] = str(self.allToOne(self.multinet.seqaln, struct.name, struct.sequenceList, int(G.nodes[i]['label'])))

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


    def exportPickle (self):

        # Creates new pickle (.pkl) file and then dumps the entire class object into the pickle file
        with open(f'{self.args.outputname}.pkl', 'wb') as pickleFile:
            pickle.dump(self, pickleFile)
