from multirin.analysis.Similarity import Similarity
from multirin.generate.MultiNetwork import MultiNetwork
from argparse import Namespace
import unittest
import numpy as np
import xarray as xr
import networkx as nx

# Define class to test the program
class testSimilarity (unittest.TestCase):

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

    def test_distanceMatrix (self):

        args = Namespace(alignmentFile='tests/data/multi_net_test/PTP-KDY.fa')
        simObject = Similarity(args)
        simObject.multinet = MultiNetwork(args)

        # Creates test array object
        simObject.multinet.array = self.multinetArray

        # Generates distance matrix
        simObject.distanceMatrix()

        # Get value of total distance in this test matrix
        val = 0
        for i in range(0,9,1):
            j = i + 9

            val = val + abs((i * 10 / 8) - (j * 10 / 17))
        
        self.assertAlmostEqual(simObject.distanceArray.loc['2SHV','1ALI'].item(), val)
        self.assertAlmostEqual(simObject.distanceArray.loc['2SHV','2SHV'].item(), 0)


if __name__ == '__main__':
    unittest.main()