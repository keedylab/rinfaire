import numpy as np
from Bio import SeqIO
import networkx as nx
from scipy.cluster.hierarchy import dendrogram, linkage, average
from scipy.spatial.distance import pdist, squareform
from matplotlib import pyplot as plt
from pyvis.network import Network
import csv

class MultiNetwork:
    
    # Class constructor
    def __init__(self,alignmentpath):
        
        self.setSeqAlignment(alignmentpath)
        self.array = np.zeros((0, self.size, self.size))
            
        #self.lookuptable = table with pdb indices and relevant information

    # Set functions

    def setSeqAlignment (self, alignmentpath):
        # Uses the SeqIO tool from Biopython to parse the fasta file
        self.seqaln = {}
        self.size = 0
        
        for record in SeqIO.parse(alignmentpath, "fasta"):
            pdbid = record.id[0:4]
            self.seqaln[pdbid] = record.seq
            
            # Makes sure that self.size is as long as the longest sequence (although they should be all the same) 
            if len(record.seq) > self.size:
                self.size = len(record.seq)
            
        # Increments by one due to how indexing starts at 0 but residue numbering starts at 1
        self.size += 1

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

    def add (self, inputAdjacencyDict, inputStructName, inputStartingResi, normalizeflag):
        
        # Creates new blank 3D array that is of size (1 x length of seq x length of seq)
        addArray = np.zeros((1, self.size, self.size))
                
        for firstResi in inputAdjacencyDict:
            for secondResi in inputAdjacencyDict[firstResi]:
                
                #print(f"Before conversion: {firstResi} {secondResi}")

                # Uses the OneToAll conversion function to map individual resi # to sequence alignment #
                updatedFirstResi = self.oneToAll(inputStructName, inputStartingResi, firstResi)
                updatedSecondResi = self.oneToAll(inputStructName, inputStartingResi, secondResi)

                #print(f"After conversion: {updatedFirstResi} {updatedSecondResi}")

                # Adds pairing to the array along with the weight
                addArray[0][updatedFirstResi][updatedSecondResi] = inputAdjacencyDict[firstResi][secondResi]['width']
                
        # if normalizeflag == True:
        #     maxvalue = 0
        #     for idx, x in np.ndenumerate(addArray):
        #         if x > maxvalue:
        #             maxvalue = x

        #     print(maxvalue)

        #     # Iterate over all elements and divide by the max (scales from 0-10)
        #     with np.nditer(addarray, op_flags=['readwrite']) as it:
        #         for x in it:
        #             x[...] = (x/maxvalue)*10
        
        # Merges new array with the main array
        self.array = np.concatenate((self.array,addArray), axis=0)
        
        # # Visualizing the added array
        # addarray2d = addarray.reshape(self.size,self.size) # Converts addarray that technically is 3D array with third dimension of size 1 back to flat 2D array
        # self.VisArray(addarray2d,addpdbname)

    def sum (self):

        sumMatrix = np.sum(self.array,axis=0)
        return sumMatrix
    
    def visualize (self, inputarray, outputname):

        # Creates graph from Numpy array
        G = nx.from_numpy_array(inputarray)
        G.remove_nodes_from(list(nx.isolates(G)))
        nx.convert_node_labels_to_integers(G)
        
        for i in G.nodes():
            G.nodes[i]['label'] = str(i)
        
        # widths = nx.get_edge_attributes(G, 'weight')

        #print("Multi Network \n", G.edges(data=True), "\n")
        
        nts = Network(notebook=True)
        
        # populates the nodes and edges data structures
        nts.from_nx(G)
        outputpath = outputname + '.html'
        nts.show(outputpath)
    

    