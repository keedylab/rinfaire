import pickle
import numpy as np
import xarray as xr
from scipy.cluster.hierarchy import dendrogram, linkage, average
from scipy.spatial.distance import pdist, squareform
import netrd

class Similarity:

    def __init__ (self, args):
        self.args = args

    def readPickle (self):
        
        # Opens pickle file
        with open(self.args.filename, 'rb') as pickleFile:
            self.multinet = pickle.load(pickleFile)

    def distanceMatrix (self):

        self.multinet.normalizeStruct()
        
        # Creates a new distance matrix in XArray
        self.distanceArray = xr.DataArray( 
            coords=dict(firstPDB=self.multinet.array.network.data, secondPDB=self.multinet.array.network.data), 
            dims=("firstPDB", "secondPDB")
        )

        inputArray = self.multinet.array
        # inputArray = self.multinet.array.where(self.multinet.array == 0, other=1)

        # Loops over every pair of PDBs in the network object
        for pdb1 in self.multinet.array.network:
            for pdb2 in self.multinet.array.network:

                # Takes the absolute value of the difference of the two slices from each PDB
                diffSlice = abs(inputArray.loc[pdb1, :, :] - inputArray.loc[pdb2, :, :])

                # Then sums the values across the array of differences at each position
                diffValue = diffSlice.sum().values

                # Appends the value of this sum to the corresponding location in the distance matrix
                self.distanceArray.loc[pdb1, pdb2] = diffValue
                #print(pdb1.values, pdb2.values, diffSum)

    def heirClustering (self):

        import matplotlib.pyplot as plt

        condensedMatrix = squareform(self.distanceArray.to_numpy())
        
        print(condensedMatrix)
        
        linkmatrix = average(condensedMatrix)
        #linkmatrix = linkage(testcondensed, method='average')
        
        print(linkmatrix)
        
        plt.figure(figsize=(50, 50))
        dn = dendrogram(linkmatrix)

        plt.tight_layout()
        filename = self.args.outputdir + 'HeirClustering'
        outputpath = f'{filename}.png'
        plt.savefig(outputpath)
    
    def visualizeMatrix (self):

        import seaborn as sns
        import matplotlib.pyplot as plt

        # Creates plot
        plt.figure(figsize=(10,10))
        sns.set(font_scale=1.5)

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