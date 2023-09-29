from multirin.analysis.SumNetwork import SumNetwork
from multirin.generate.MultiNetwork import MultiNetwork
from argparse import Namespace
import unittest
import numpy as np
import xarray as xr

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

    def test_calculateSum (self):

        args = Namespace(alignmentFile='tests/data/multi_net_test/PTP-KDY.fa', no_scale_sum_network=True, remove_weak_edges=None)
        sumNetwork = SumNetwork(args)
        sumNetwork.multinet = MultiNetwork(args)

        # Creates test array object
        sumNetwork.multinet.array = self.multinetArray

        # Creates sum array and then scales it
        sumNetwork.calculateSum() 

        # Assert that values are equal
        self.assertEqual(sumNetwork.sumArray[0][0].item(), (sumNetwork.multinet.array[0][0][0].item() + sumNetwork.multinet.array[1][0][0].item()))
        self.assertEqual(sumNetwork.sumArray[2][2].item(), (sumNetwork.multinet.array[0][2][2].item() + sumNetwork.multinet.array[1][2][2].item()))

    def test_scaleSumNetwork (self):

        args = Namespace(alignmentFile='tests/data/multi_net_test/PTP-KDY.fa', no_scale_sum_network=False, sum_network_scale=20, remove_weak_edges=None)
        sumNetwork = SumNetwork(args)
        sumNetwork.multinet = MultiNetwork(args)

        # Creates test array object
        sumNetwork.multinet.array = self.multinetArray

        # Creates sum array and then scales it
        sumNetwork.calculateSum() 

        # Assert that values are equal
        maxValue = 25 # This is the maximum value in the array before scaling
        self.assertEqual(sumNetwork.sumArray[0][0].item(), ((sumNetwork.multinet.array[0][0][0].item() + sumNetwork.multinet.array[1][0][0].item()) / maxValue) * sumNetwork.args.sum_network_scale)
        self.assertEqual(sumNetwork.sumArray[2][2].item(), ((sumNetwork.multinet.array[0][2][2].item() + sumNetwork.multinet.array[1][2][2].item()) / maxValue) * sumNetwork.args.sum_network_scale)

    def test_removeWeakEdges (self):

        args = Namespace(alignmentFile='tests/data/multi_net_test/PTP-KDY.fa', no_scale_sum_network=False, sum_network_scale=20, remove_weak_edges=20)
        sumNetwork = SumNetwork(args)
        sumNetwork.multinet = MultiNetwork(args)

        # Creates test array object
        sumNetwork.multinet.array = self.multinetArray

        # Creates sum array and then scales it
        # This also will remove weak edges below 20% cutoff since the flag is on
        sumNetwork.calculateSum()

        # Original array's maximum value (after scaling it should be 20)
        maxValue = 20

        # Gets array of all values that are larger than zero but less than the cutoff
        # Ideally this should be empty
        cutoffArray = sumNetwork.sumArray.where((sumNetwork.sumArray < 0.2 * maxValue) & (sumNetwork.sumArray > 0), drop=True)

        # Asserts that this array is empty (since there should be no values after removal that fit criteria)
        self.assertEqual(len(cutoffArray), 0)
 
if __name__ == '__main__':
    unittest.main()