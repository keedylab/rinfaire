from multirin.analysis.Covariance import Covariance
from multirin.generate.MultiNetwork import MultiNetwork
from argparse import Namespace
import unittest
import numpy as np
import xarray as xr
import networkx as nx

# Define class to test the program
class testSumNetwork (unittest.TestCase):

    @classmethod
    def setUpClass (self):

        # Sets up test multidimensional array that is common for all tests
        structList = ["2SHV","1ALI"]
        self.multinetArray = xr.DataArray(
            np.arange(18).reshape(2, 3, 3), 
            coords=dict(network=structList, firstResi=range(3), secondResi=range(3)), 
            dims=("network", "firstResi", "secondResi"))
        
        structList2 = ['a','b','c','d','e','f','g','h','i','j']
        self.multinetArray2 = xr.DataArray(
            np.random.rand(10, 3, 3), 
            coords=dict(network=structList2, firstResi=range(3), secondResi=range(3)), 
            dims=("network", "firstResi", "secondResi"))

    def test_flatten (self):

        args = Namespace(alignmentFile='tests/data/multi_net_test/PTP-KDY.fa', no_scale_sum_network=True, remove_weak_edges=None)
        covObject = Covariance(args)
        covObject.multinet = MultiNetwork(args)

        # Creates test array object
        covObject.multinet.array = self.multinetArray

        # Flattens the covariance array
        flattenedArray = covObject.flatten()

        # Assert that values are equal
        self.assertEqual(flattenedArray.sizes['network'], 2)
        self.assertEqual(flattenedArray.sizes['resiPair'], 9)

    def test_calculateCovarianceByResiPair (self):

        args = Namespace(alignmentFile='tests/data/multi_net_test/PTP-KDY.fa', no_scale_sum_network=True, remove_weak_edges=None)
        covObject = Covariance(args)
        covObject.multinet = MultiNetwork(args)

        # Creates test array object
        covObject.multinet.array = self.multinetArray

        # Calculates the covariance by resi pair
        covObject.calculateCovarianceByResiPair()
        
        # Test case of two different vectors
        # ResiPair Vectors are: (0,0) and (1,0)
        # (0,0) Values: [0,9]
        # (1,0) Values: [3,12]
        meanX1 = (9+0)/2
        meanY1 = (3+12)/2
        covariance1 = ((0-meanX1)*(3-meanY1) + (9-meanX1)*(12-meanY1)) / (2-1)
        
        # Test case on the diagonal where the covariance = variance
        # ResiPair Vectors are: (1,2) and (1,2)
        # (1,2) Values: [5,14]
        meanX2 = (5+14)/2
        covariance2 = ((5-meanX2)*(5-meanX2) + (14-meanX2)*(14-meanX2)) / (2-1)

        self.assertEqual(covObject.covarianceArray.loc[[(0,0)], [(1,0)]].item(), covariance1)
        self.assertEqual(covObject.covarianceArray.loc[[(1,2)], [(1,2)]].item(), covariance2)

    def test_calculateCorrelationByResiPair (self):

        args = Namespace(alignmentFile='tests/data/multi_net_test/PTP-KDY.fa', no_scale_sum_network=True, remove_weak_edges=None)
        covObject = Covariance(args)
        covObject.multinet = MultiNetwork(args)

        # Creates test array object
        covObject.multinet.array = self.multinetArray2

        # Calculates the covariance by resi pair
        covObject.calculateCorrelationByResiPair()

        # Confirms that the matrix is symmetric (transpose = original)
        self.assertTrue((covObject.correlationArray.T == covObject.correlationArray).all())

        # Confirms diagonal elements = 1
        diagValues = np.diag(covObject.correlationArray)
        diagValues = np.rint(diagValues).astype(int)
        self.assertTrue((diagValues == 1).all())

if __name__ == '__main__':
    unittest.main()