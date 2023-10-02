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

    def test_flatten (self):

        args = Namespace(alignmentFile='tests/data/multi_net_test/PTP-KDY.fa', no_scale_sum_network=True, remove_weak_edges=None)
        covObject = Covariance(args)
        covObject.multinet = MultiNetwork(args)

        # Creates test array object
        covObject.multinet.array = self.multinetArray

        # Creates sum array and then scales it
        covObject.flatten()

        # Assert that values are equal
        self.assertEqual(covObject.flattenMulti.sizes['network'], 2)
        self.assertEqual(covObject.flattenMulti.sizes['resiPair'], 9)

if __name__ == '__main__':
    unittest.main()