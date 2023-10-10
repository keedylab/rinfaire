import xarray as xr
import numpy as np
import networkx as nx
from pyvis.network import Network
import logging
import pickle
import seaborn as sns
import matplotlib.pyplot as plt

class Covariance:

    def __init__ (self, args):
        self.args = args

    def readPickle (self):
        
        # Opens pickle file
        with open(self.args.filename, 'rb') as pickleFile:
            self.multinet = pickle.load(pickleFile)

    def flatten (self):
        return self.multinet.array.stack(resiPair=("firstResi", "secondResi"))

    def calculateCovarianceByResiPair (self):
        
        # Flattens the array so that residues are now represented on a single axis as pairs such as (resi i, j)
        flattenedArray = self.flatten()

        # Converts the flattened array into a numpy array
        # Then uses numpy to calculate the sample covariance
        # - rowvar = false means that columns (residue pairs) are the variables instead of the rows (network/PDB)
        covarianceArrayNP = np.cov(flattenedArray.to_numpy(), rowvar=False)

        # Then wraps the numpy array back into an XArray dataset with the same labels as before
        # Now a square matrix with dimension of: number of residue pairs x number of residue pairs
        self.covarianceArray = xr.DataArray(
            covarianceArrayNP, 
            coords=dict(firstPair=flattenedArray.resiPair.data, secondPair=flattenedArray.resiPair.data), 
            dims=("firstPair", "secondPair")
        )

    def visualize (self):
        
        plt.figure(figsize=(10,10))
        sns.set(font_scale=1.5)
        hm = sns.heatmap(self.covarianceArray)
            # cbar=True,
            # square=True,
            # fmt='.2f',
            # annot_kws={'size': 12}
            # yticklabels=cols,
            # xticklabels=cols
        plt.title('Covariance matrix')
        plt.tight_layout()
        
        outputpath = f'{self.args.outputname}.png'
        plt.savefig(outputpath)

    def exportPickle (self):

        # Creates new pickle (.pkl) file and then dumps the entire class object into the pickle file
        with open(f'{self.args.outputname}.pkl', 'wb') as pickleFile:
            pickle.dump(self, pickleFile)
