import xarray as xr
import numpy as np
import networkx as nx
from pyvis.network import Network
import logging
import pickle
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import squareform
from sklearn.preprocessing import StandardScaler

class Covariance:

    def __init__ (self, args):
        self.args = args

    def readPickle (self):
        
        # Opens pickle file
        with open(self.args.filename, 'rb') as pickleFile:
            self.multinet = pickle.load(pickleFile)

    def removeWeakEdges (self, sumArray):
        
        # Sorts the array and finds the maximum index to keep (by taking size of array * percent to cutoff)
        # Then gets value at this index
        sortedArray = np.sort(sumArray, axis=None)
        sortedArray = sortedArray[sortedArray != 0]
        maxIndex = round(sortedArray.size * (self.args.remove_weak_edges / 100))
        cutoffValue = sortedArray[maxIndex]

        # Finds entries in the array where they are less than the cutoff threshold
        # Then it replaces them with 0.0
        # If not, then it keeps the original value
        sumArray = xr.where(sumArray < cutoffValue, 0.0, sumArray)
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

    def calculateCovarianceByResiPair (self, scaleFlag):
        
        # Flattens the array so that residues are now represented on a single axis as pairs such as (resi i, j)
        flattenedArray = self.flatten()

        # Does standard scaling for PCA implementation
        if scaleFlag == True:
            flattenedArray = (flattenedArray - flattenedArray.mean(axis = 0)) / flattenedArray.std(axis = 0)

        # Converts the flattened array into a numpy array
        # Then uses numpy to calculate the sample covariance
        # rowvar = false means that columns (residue pairs) are the variables instead of the rows (network/PDB)
        self.covarianceArrayNP = np.cov(flattenedArray.to_numpy(), rowvar=False)

        # Then wraps the numpy array back into an XArray dataset with the same labels as before
        # Now a square matrix with dimension of: number of residue pairs x number of residue pairs
        self.covarianceArray = xr.DataArray(
            self.covarianceArrayNP, 
            coords=dict(firstPair=flattenedArray.resiPair.data, secondPair=flattenedArray.resiPair.data), 
            dims=("firstPair", "secondPair")
        )

        # Gets largest and smallest covariance values that are not on the diagonal (variance with itself)
        # Solution from: https://stackoverflow.com/questions/29394377/minimum-of-numpy-array-ignoring-diagonal
        covarianceArrayNP_Print = self.covarianceArrayNP.copy()
        np.fill_diagonal(covarianceArrayNP_Print, -np.inf)
        max_value = covarianceArrayNP_Print.max()
        np.fill_diagonal(covarianceArrayNP_Print, np.inf)
        min_value = covarianceArrayNP_Print.min()

        print(f"Minimum covariance (w/o diagonal entries): {str(min_value)}")
        print(f"Maximum covariance (w/o diagonal entries): {str(max_value)}")

    def calculateCorrelationByResiPair (self):
        
        # Flattens the array so that residues are now represented on a single axis as pairs such as (resi i, j)
        flattenedArray = self.flatten()

        # Converts the flattened array into a numpy array
        # Then uses numpy to calculate the Pearson correlation coefficient
        # rowvar = false means that columns (residue pairs) are the variables instead of the rows (network/PDB)
        self.correlationArrayNP = np.corrcoef(flattenedArray.to_numpy(), rowvar=False)
        self.correlationArrayNP = np.around(self.correlationArrayNP, decimals=3) # Rounds this to nearest thousandth

        # Then wraps the numpy array back into an XArray dataset with the same labels as before
        # Now a square matrix with dimension of: number of residue pairs x number of residue pairs
        self.correlationArray = xr.DataArray(
            self.correlationArrayNP, 
            coords=dict(firstPair=flattenedArray.resiPair.data, secondPair=flattenedArray.resiPair.data), 
            dims=("firstPair", "secondPair")
        )

        # Gets largest and smallest correlation coefficients that are not on the diagonal (correlation with itself)
        # Solution from: https://stackoverflow.com/questions/29394377/minimum-of-numpy-array-ignoring-diagonal
        correlationArrayNP_Print = self.correlationArrayNP.copy()
        np.fill_diagonal(correlationArrayNP_Print, -np.inf)
        max_value = correlationArrayNP_Print.max()
        np.fill_diagonal(correlationArrayNP_Print, np.inf)
        min_value = correlationArrayNP_Print.min()

        print(f"Minimum Pearson correlation coefficient (w/o diagonal entries): {str(min_value)}")
        print(f"Maximum Pearson correlation coefficient (w/o diagonal entries): {str(max_value)}")

    def runPCA (self):

        # Gets eigenvalues and eigenvectors of the covariance matrix
        # This allows us to diagonalize the covariance array
        # Uses np.linalg.eigh because we know it's a symmetric matrix and therefore has real eigenvalues
        eigenvalues, eigenvectors = np.linalg.eigh(self.covarianceArrayNP)

        # np.argsort can only provide lowest to highest; use [::-1] to reverse the list
        order_of_importance = np.argsort(eigenvalues)[::-1] 

        # utilize the sort order to sort eigenvalues and eigenvectors
        sorted_eigenvalues = eigenvalues[order_of_importance]
        sorted_eigenvectors = eigenvectors[:,order_of_importance] # sort the columns

        # use sorted_eigenvalues to ensure the explained variances correspond to the eigenvectors
        explained_variance = sorted_eigenvalues / np.sum(sorted_eigenvalues)
        
        # Gets the cumulative sum of % explained variance for increasing numbers of principal components (PCs)
        cumSumVariance = np.cumsum(explained_variance)

        # Gets the number of PCs required to explain 95% of the total variance
        PCIndex = 0
        for PC in cumSumVariance:
            
            PCIndex += 1
            if PC >= 0.95:
                break

        if self.args.output_variance_plot == True:
            # Plots the cumulative sum for increasing numbers of PCs
            plt.plot(cumSumVariance)
            filename = self.args.outputdir + 'PCA_ExplainedVariance'
            outputpath = f'{filename}.png'
            plt.savefig(outputpath)

        # Multiply original matrix by eigenvector matrix
        # This projects each datapoint onto eigenvectors that are chosen
        covarianceArrayPCA = np.matmul(self.covarianceArrayNP, sorted_eigenvectors[:,:PCIndex])
        print(covarianceArrayPCA.max(), covarianceArrayPCA.min())

        return covarianceArrayPCA

    def clusterCorrMatrix (self):

        # Overall solution from: https://www.kaggle.com/code/sgalella/correlation-heatmaps-with-hierarchical-clustering
        # First creates a distance/dissimilarity matrix by taking 1 minus the abs val of each element in the correlation array
        # This makes it so that highly correlated Pearson coefficient values (1 or -1) have a distance of 0
        # And that Pearson coefficient values with no correlation (0) have a distance of 1 (maximum distance in this case)
        dissimilarity = 1 - abs(self.correlationArrayNP)
        print(dissimilarity)

        # Then it creates a linkage matrix
        linkageMatrix = linkage(squareform(dissimilarity), 'average')

        plt.figure(figsize=(12,5))
        dendrogramCorr = dendrogram(linkageMatrix, labels=self.correlationArray['firstPair'], orientation='top', 
                leaf_rotation=90)
        
        filename = self.args.outputdir + 'CorrDendrogram'
        outputpath = f'{filename}.png'
        plt.savefig(outputpath)

    def visualizeMatrix (self, covOrCorrFlag):

        # Creates plot
        plt.figure(figsize=(10,10))
        sns.set(font_scale=1.5)

        # Chooses the input array depending on if it's the covariance or correlation matrix
        if covOrCorrFlag == 'covariance':
            visArray = self.covarianceArray
            plt.title('Covariance matrix')
            filename = self.args.outputdir + 'Covariance'
        else:
            visArray = self.correlationArray
            plt.title('Correlation matrix')
            filename = self.args.outputdir + 'Correlation'
        
        # Creates a heatmap for the matrix
        hm = sns.heatmap(visArray)
            # cbar=True,
            # square=True,
            # fmt='.2f',
            # annot_kws={'size': 12}
            # yticklabels=cols,
            # xticklabels=cols
        
        # Plots the heatmap and saves it to the specified directory
        plt.tight_layout()
        outputpath = f'{filename}.png'
        plt.savefig(outputpath)

        # plt.figure(figsize=(10,10))
        # hmCluster = sns.clustermap(visArray)
        # outputpath = f'{filename}_Cluster.png'
        # plt.savefig(outputpath)

    def exportPickle (self):

        # Creates new pickle (.pkl) file and then dumps the entire class object into the pickle file
        with open(f'{self.args.outputdir + "Covariance"}.pkl', 'wb') as pickleFile:
            pickle.dump(self, pickleFile)
