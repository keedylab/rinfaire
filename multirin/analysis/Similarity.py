import pickle
import numpy as np
import xarray as xr
from scipy.cluster.hierarchy import dendrogram, linkage, average
from scipy.spatial.distance import pdist, squareform

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

        # Loops over every pair of PDBs in the network object
        for pdb1 in self.multinet.array.network:
            for pdb2 in self.multinet.array.network:

                # Takes the absolute value of the difference of the two slices from each PDB
                diffSlice = abs(self.multinet.array.loc[pdb1, :, :] - self.multinet.array.loc[pdb2, :, :])

                # Then sums the values across the array of differences at each position
                diffSum = diffSlice.sum().values

                # Appends the value of this sum to the corresponding location in the distance matrix
                self.distanceArray.loc[pdb1, pdb2] = diffSum
                #print(pdb1.values, pdb2.values, diffSum)
    
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