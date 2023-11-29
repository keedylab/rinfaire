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

    def removeWeakEdges (self, sumArray):
        
        # Get maximum value across the entire array
        # Then multiply that by the percent cutoff threshold specified by user
        percentCutoffValue = sumArray.max() * (self.args.remove_weak_edges / 100)

        # Finds entries in the array where they are less than the percent cutoff threshold
        # Then it replaces them with 0.0
        # If not, then it keeps the original value
        sumArray = xr.where(sumArray < percentCutoffValue, 0.0, sumArray)
        return sumArray

    def flatten (self):

        # First stacks the array by creating a combined coordinate resiPair that combines firstResi and secondResi
        stackedArray = self.multinet.array.stack(resiPair=("firstResi", "secondResi"))
        stackedArray = stackedArray.dropna(dim="network", how="any")
        
        # Then calculates the sum of the stackedArray across the network dimension (ie the sum for each resiPair)
        sumArray = stackedArray.sum(dim="network")

        # Removes weak edges if option is specified
        # Makes entries in array 0 if they fall below the threshold -> entries that are = 0 get removed in next step
        if self.args.remove_weak_edges != None:
            sumArray = self.removeWeakEdges(sumArray)

        # Then drops indices where the corresponding index in the sumArray is = 0 (meaning that there are no edges for that pair)
        # Creates a stackedArray with only columns with resiPairs that exist in the network
        stackedArray = stackedArray.where(sumArray != 0, drop=True)
        return stackedArray

    def calculateCovarianceByResiPair (self):
        
        # Flattens the array so that residues are now represented on a single axis as pairs such as (resi i, j)
        flattenedArray = self.flatten()

        # Converts the flattened array into a numpy array
        # Then uses numpy to calculate the sample covariance
        # rowvar = false means that columns (residue pairs) are the variables instead of the rows (network/PDB)
        covarianceArrayNP = np.cov(flattenedArray.to_numpy(), rowvar=False)

        # Then wraps the numpy array back into an XArray dataset with the same labels as before
        # Now a square matrix with dimension of: number of residue pairs x number of residue pairs
        self.covarianceArray = xr.DataArray(
            covarianceArrayNP, 
            coords=dict(firstPair=flattenedArray.resiPair.data, secondPair=flattenedArray.resiPair.data), 
            dims=("firstPair", "secondPair")
        )

        print(self.covarianceArray)

    def calculateCorrelationByResiPair (self):
        
        # Flattens the array so that residues are now represented on a single axis as pairs such as (resi i, j)
        flattenedArray = self.flatten()

        # Converts the flattened array into a numpy array
        # Then uses numpy to calculate the Pearson correlation coefficient
        # rowvar = false means that columns (residue pairs) are the variables instead of the rows (network/PDB)
        correlationArrayNP = np.corrcoef(flattenedArray.to_numpy(), rowvar=False)

        # Then wraps the numpy array back into an XArray dataset with the same labels as before
        # Now a square matrix with dimension of: number of residue pairs x number of residue pairs
        self.correlationArray = xr.DataArray(
            correlationArrayNP, 
            coords=dict(firstPair=flattenedArray.resiPair.data, secondPair=flattenedArray.resiPair.data), 
            dims=("firstPair", "secondPair")
        )

    def visualize (self):
        
        plt.figure(figsize=(10,10))
        sns.set(font_scale=1.5)
        hm = sns.heatmap(self.correlationArray)
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
