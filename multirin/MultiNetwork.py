# import numpy as np
import xarray as xr
from Bio import SeqIO
import networkx as nx
from scipy.cluster.hierarchy import dendrogram, linkage, average
from scipy.spatial.distance import pdist, squareform
from matplotlib import pyplot as plt
from pyvis.network import Network
import logging

class MultiNetwork:
    
    # Class constructor
    def __init__ (self, args):
        
        self.setSeqAlignment(args.alignmentFile)
        self.array = None
        # self.array = np.zeros((0, self.size, self.size))
        self.args = args

        # if args.labels is not None:
        #     self.setMetaData(args.labels)
            
        #self.lookuptable = table with pdb indices and relevant information

    # Set functions

    def setSeqAlignment (self, alignmentPath):
        # Uses the SeqIO tool from Biopython to parse the fasta file
        self.seqaln = {}
        self.size = 0
        
        for record in SeqIO.parse(alignmentPath, "fasta"):
            pdbid = record.id[0:4]
            self.seqaln[pdbid] = record.seq
            
            # Makes sure that self.size is as long as the longest sequence (although they should be all the same) 
            if len(record.seq) > self.size:
                self.size = len(record.seq)
            
        # Increments by one due to how indexing starts at 0 but residue numbering starts at 1
        self.size += 1

    # TODO: Create metadata function
    # def setMetaData (self, csvFile):

    #     self.metadata = None

    def oneToAll (self, seqID, startingResi, seqResidue):
        
        mainCount = 0
        seqCount = startingResi - 1 # Subtracts 1 b/c it has not iterated over that residue yet in its traversal of the sequence
   
        #if self.seqaln[seq_id][seq_residue] == '-':
        #    return None;
        
        # Iterates over each residue in the sequence of the structure being queried
        for i in self.seqaln[seqID]:

            # Increments the main count that counts the total number of positions it has passed on the alignment
            mainCount = mainCount + 1

            # If there is a residue present (not just a - used as a placeholder)
            # Then increment the seq count that counts the total number of actual residue positions that have been passed
            if i != '-':
                seqCount = seqCount + 1

            # Condition when it reaches the residue in question 
            # (since seq count represents the position in the individual sequence and not the full alignment)
            if seqCount == int(seqResidue) and i != '-':
                
                #print("Final Main count:", mainCount," Seq count: ", seqCount, 'Resi:', i)

                # Returns the main count which is the position on the full alignment
                return(mainCount)

            #print("Main count:", mainCount," Seq count: ", seqCount, 'Resi:', i)

    def add (self, inputAdjacencyDict, inputStructName, inputStartingResi):

        # Loops through all pairings                
        for firstResi in inputAdjacencyDict:
            for secondResi in inputAdjacencyDict[firstResi]:
                
                #print(f"Before conversion: {firstResi} {secondResi}")

                # Uses the OneToAll conversion function to map individual resi # to sequence alignment #
                updatedFirstResi = self.oneToAll(inputStructName[0:4], inputStartingResi, firstResi)
                updatedSecondResi = self.oneToAll(inputStructName[0:4], inputStartingResi, secondResi)

                #print(f"After conversion: {updatedFirstResi} {updatedSecondResi}")

                # Adds pairing to the array along with the weight
                self.array.loc[inputStructName, updatedFirstResi, updatedSecondResi] = inputAdjacencyDict[firstResi][secondResi]['weight']

    # TODO: Update unit test to make sure this function works
    def normalizeStruct (self):

        # Creates a vector of maximum values across the first and second residue
        # Essentially a maximum value for each network
        maxValues = self.array.max(dim=['firstResi','secondResi'])

        # Then divides each network by the corresponding value in the vector of max values
        # Scales from 0 - 10
        self.array = (self.array / maxValues[:]) * 10
    
    # TODO: Update unit test to make sure this function works
    def addNetworks (self, networkList):

        # Creates list to represent all the structures
        structList = []
        for net in networkList:
            structList.append(net.struct.getName())

        print(structList)

        # Creates new blank 3D array that is of size (number of structures x length of seq x length of seq)
        self.array = xr.DataArray(
            0, 
            coords=dict(network=structList, firstResi=range(self.size), secondResi=range(self.size)), 
            dims=("network", "firstResi", "secondResi")
        )

        # Iterates over the network list and adds values from each adjacency matrix
        for net in networkList:
            self.add(
                net.convertToAdjacency(), 
                net.struct.getName(), 
                net.struct.getFirstResi()
            )

        # Normalization of edge weights relative to the whole structure being added
        # All values are scaled from (0) to the max edge weight present (10)
        if self.args.no_norm_struct == False:
            self.normalizeStruct()

        print(self.array)
        
    # TODO: Update unit test to make sure this function works
    def sum (self):

        sumMatrix = self.array.sum(dim="network")
        return sumMatrix
    
    def visualize (self, inputarray, outputname):

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
        outputpath = outputname + '.html'
        nts.show(outputpath)
    

    