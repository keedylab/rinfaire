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

    def test_removeSubGraphs (self):

        args = Namespace(remove_subgraphs=4)
        sumNetwork = SumNetwork(args)

        # Creates empty networkX graph and then populates it with edges of different weights
        graph = nx.Graph()
        graph.add_edge(2, 1, weight=4)
        graph.add_edge(2, 3, weight=3)
        graph.add_edge(2, 5, weight=2)
        graph.add_edge(7, 9, weight=2)
        graph.add_edge(9, 10, weight=2)
        graph.add_edge(11, 12, weight=5)

        # Removes subgraphs
        graphRemoved = sumNetwork.removeSubGraphs(graph)

        # Checks what nodes are present
        self.assertEqual(sorted(list(graphRemoved.nodes)), sorted([1,2,3,5]))

    def test_seqToRef (self):

        args = Namespace(alignmentFile='tests/data/multi_net_test/PTP-KDY.fa', no_scale_sum_network=False, sum_network_scale=20, remove_weak_edges=None, seq_to_ref='tests/data/sum_net_test/6B8Z_qFit_chainA.pdb')
        sumNetwork = SumNetwork(args)
        sumNetwork.multinet = MultiNetwork(args)

        # Creates test array object
        sumNetwork.multinet.array = self.multinetArray

        # Creates sum array and then scales it
        sumNetwork.seqToRef() 

        # Assert that values are equal
        # maxValue = 25 # This is the maximum value in the array before scaling
        # self.assertEqual(sumNetwork.sumArray[0][0].item(), ((sumNetwork.multinet.array[0][0][0].item() + sumNetwork.multinet.array[1][0][0].item()) / maxValue) * sumNetwork.args.sum_network_scale)
        # self.assertEqual(sumNetwork.sumArray[2][2].item(), ((sumNetwork.multinet.array[0][2][2].item() + sumNetwork.multinet.array[1][2][2].item()) / maxValue) * sumNetwork.args.sum_network_scale)

        

if __name__ == '__main__':
    unittest.main()