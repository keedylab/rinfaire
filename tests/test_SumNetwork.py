from multirin.analysis.SumNetwork import SumNetwork
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

        args = Namespace(alignmentFile='tests/data/multi_net_test/PTP-KDY.fa', no_scale_sum_network=False, sum_network_scale=20, remove_weak_edges=None)
        sumNetwork = SumNetwork(args)
        sumNetwork.multinet = MultiNetwork(args)

        # Creates test array object
        sumNetwork.multinet.array = self.multinetArray
        sumNetwork.calculateSum()

        sumNonRemoved = sumNetwork.sumArray.to_numpy().flatten()
        lenNonRemoved = len(sumNonRemoved[sumNonRemoved != 0])

        # Creates another test object for the removed edges case
        args = Namespace(alignmentFile='tests/data/multi_net_test/PTP-KDY.fa', no_scale_sum_network=False, sum_network_scale=20, remove_weak_edges=20)
        sumNetworkRemoved = SumNetwork(args)
        sumNetworkRemoved.multinet = MultiNetwork(args)

        # Creates test array object
        sumNetworkRemoved.multinet.array = self.multinetArray
        sumNetworkRemoved.calculateSum()

        sumRemoved = sumNetworkRemoved.sumArray.to_numpy().flatten()
        lenRemoved = len(sumRemoved[sumRemoved != 0])

        self.assertEqual(lenRemoved, round(lenNonRemoved * ((100 - args.remove_weak_edges) / 100)))

    def test_resizeByDegree (self):

        args = Namespace(resize_by_degree_scale=2)
        sumNetwork = SumNetwork(args)

        # Creates empty networkX graph and then populates it with edges of different weights
        graph = nx.Graph()
        graph.add_edge(2, 1, weight=4.7)
        graph.add_edge(2, 3, weight=3.3)
        graph.add_edge(2, 5, weight=2.0)

        # Resizes the graph
        graphResized = sumNetwork.resizeByDegree(graph)

        # Checks if the node sizes are what they should be
        self.assertEqual(graphResized.nodes[2]['size'], 20)
        self.assertEqual(graphResized.nodes[1]['size'], 9.4)
        self.assertEqual(graphResized.nodes[3]['size'], 6.6)
        self.assertEqual(graphResized.nodes[5]['size'], 4.0)
 
if __name__ == '__main__':
    unittest.main()