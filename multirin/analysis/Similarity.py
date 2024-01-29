import pickle
import numpy as np
import xarray as xr
from scipy.cluster.hierarchy import dendrogram, linkage, average, fcluster
from scipy.spatial.distance import pdist, squareform
import netrd
import networkx as nx

class Similarity:

    def __init__ (self, args):
        self.args = args

    def readPickle (self):
        
        # Opens pickle file
        with open(self.args.filename, 'rb') as pickleFile:
            self.multinet = pickle.load(pickleFile)

    def distanceMatrix (self):

        self.multinet.normalizeStruct()
        
        inputArray = self.multinet.array
        inputArray = inputArray.fillna(0)

        # Creates dictionary of networks from each adjacency matrix in the input array
        self.networkDict = {}
        
        for pdb in self.multinet.array.network:
            
            pdbName = pdb.network.data.tolist()
            netAdd = nx.from_numpy_array(inputArray.loc[pdb, :, :].to_numpy())

            # Removes networks with fewer than 5 edges for this analysis
            if netAdd.number_of_edges() >= 5:
                self.networkDict[pdbName] = netAdd

        # Creates a new distance matrix in XArray
        networkDictList = list(self.networkDict.keys())

        self.distanceArray = xr.DataArray( 
            coords=dict(firstPDB = networkDictList, secondPDB = networkDictList), 
            dims=("firstPDB", "secondPDB")
        )
        
        # Loops over every pair of PDBs in the network object
        searchedPDBs = []

        for pdb1 in networkDictList:
            for pdb2 in networkDictList:

                # Sorts pair of PDBs so we can check if we searched it or not later on
                sortedPair = sorted((pdb1, pdb2))

                # Sets the diagonal values of the matrix to be zero (distance to self = 0)
                if (pdb1 == pdb2):
                    self.distanceArray.loc[pdb1, pdb2] = 0

                # Asserts that the two PDBs are not the same and that we haven't searched this pair already
                elif (pdb1 != pdb2) and (sortedPair not in searchedPDBs):

                    # Chooses the distance metric based on the flag provided
                    if self.args.distance_metric == 'jaccard':
                        distanceObject = netrd.distance.JaccardDistance()
                    elif self.args.distance_metric == 'deltacon':
                        distanceObject = netrd.distance.DeltaCon()
                        
                    # Computes the distance between these two networks
                    graphDistance = distanceObject.dist(self.networkDict[pdb1], self.networkDict[pdb2])
                    print(pdb1, pdb2, graphDistance)

                    # Add tuple to list of searched PDB pairs
                    searchedPDBs.append(sortedPair)

                    # Appends the value of this sum to the corresponding location in the distance matrix
                    self.distanceArray.loc[pdb1, pdb2] = graphDistance
                    self.distanceArray.loc[pdb2, pdb1] = graphDistance

    def heirClustering (self):

        networkDictList = list(self.networkDict.keys())

        condensedMatrix = squareform(self.distanceArray.to_numpy())        
        linkmatrix = average(condensedMatrix)
        #linkmatrix = linkage(testcondensed, method='average')

        labels = fcluster(linkmatrix, 0.7, criterion='distance')
        clusters = self.classifyClusters(labels, networkDictList)

        # Removes clusters with only one element
        listToDel = []
        for key in clusters.keys():
            if len(clusters[key]) == 1:
                listToDel.append(key)
        
        for key in listToDel:
            del clusters[key]

        print(clusters)

        # Keep the indices to sort labels
        labelsOrder = np.argsort(labels)

        # Then reorders correlation matrix
        self.distanceArray = self.distanceArray[labelsOrder, :][:, labelsOrder]

        # Visualizes the dendrogram from heirarchical clustering
        self.visualizeDendrogram(linkmatrix)

    def visualizeDendrogram (self,linkageMatrix):

        import matplotlib.pyplot as plt

        networkDictList = list(self.networkDict.keys())
        
        plt.figure(figsize=(15, 15))
        dn = dendrogram(linkageMatrix, labels = networkDictList)

        plt.tight_layout()
        filename = self.args.outputdir + 'HeirClustering'
        outputpath = f'{filename}.png'
        plt.savefig(outputpath)

    def classifyClusters (self, labels, data):

        # Creates a dictionary of lists with all the edges in each cluster
        clusterDict = {}
        
        # Loops through the array of labels that correspond to each entry resi pair entry
        idxCounter = 0
        for labelValue in labels:

            # If it's a new label value then append a new list to the dictionary
            if labelValue not in clusterDict.keys():
                clusterDict[labelValue] = []
            
            # If not then just append the resi pair to the existing list for that dictionary key
            clusterDict[labelValue].append(data[idxCounter])

            idxCounter += 1

        return clusterDict
    
    def visualizeMatrix (self):

        import seaborn as sns
        import matplotlib.pyplot as plt

        # Creates plot
        plt.figure(figsize=(15,15))
        sns.set(font_scale=1.0)

        visArray = self.distanceArray
        plt.title('Distance matrix')
        filename = self.args.outputdir + 'DistanceMatrix'
        
        # Creates a heatmap for the matrix
        hm = sns.heatmap(visArray,
            cbar=True,
            square=True,
            # fmt='.2f',
            # annot_kws={'size': 12}
            # yticklabels=cols,
            # xticklabels=cols
        )

        # Plots the heatmap and saves it to the specified directory
        plt.tight_layout()
        outputpath = f'{filename}.png'
        plt.savefig(outputpath)