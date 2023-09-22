from multirin import MultiNetwork 
from argparse import Namespace
import unittest
import numpy as np
import xarray as xr

# Define class to test the program
class testGenerateMultiNetwork (unittest.TestCase):

    @classmethod
    def setUpClass (self):

        # Creates three different dict of dicts for testing connections
        self.Dict_2SHV = {
            1: {3: {'weight': 0.2}, 5: {'weight': 0.1}},
            3: {1: {'weight': 0.2}},
            5: {1: {'weight': 0.1}}
        }

        self.Dict_1ALI = {
            5: {7: {'weight': 0.2}, 8: {'weight': 0.1}},
            7: {5: {'weight': 0.2}},
            8: {5: {'weight': 0.1}}
        }

        self.Dict_2SJR = {
            211: {213: {'weight': 0.2}, 215: {'weight': 0.1}},
            213: {211: {'weight': 0.2}},
            215: {211: {'weight': 0.1}},
            220: {221: {'weight': 0.4}},
            221: {220: {'weight': 0.4}}
        }

    # Function to test residue conversion function
    def test_oneToAll (self):

        args = Namespace(alignmentFile='tests/data/multi_net_test/PTP-KDY.fa', no_norm_struct=False)
        multi = MultiNetwork(args)

        seqList1 = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
        seqList2 = [0, 1, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]
        seqList3 = [0, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222]

        # Tests with inputs of:
        #   Fake PDB entry: 2SHV, 3LJR
        #   Starting residue: 5, 200
        #   Residue to convert: 10, 210
        result1 = multi.oneToAll('1ALI', seqList1, 5)
        result2 = multi.oneToAll('2SHV', seqList2, 8)
        result3 = multi.oneToAll('3LJR', seqList3, 215)
        self.assertEqual(result1, 7)
        self.assertEqual(result2, 7)
        self.assertEqual(result3, 7)

    def test_add_withoutNorm (self):

        args = Namespace(alignmentFile='tests/data/multi_net_test/PTP-KDY.fa', no_norm_struct=True)
        multi = MultiNetwork(args)
        structList = ["2SHV","1ALI","2SJR"]

        multi.array = xr.DataArray(
            0.0, 
            coords=dict(network=structList, firstResi=range(multi.size), secondResi=range(multi.size)), 
            dims=("network", "firstResi", "secondResi")
        )

        seqList1 = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
        seqList2 = [0, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]
        seqList3 = [0, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222]

        multi.add(self.Dict_2SHV, "2SHV", seqList1)
        self.assertEqual(multi.array.loc["2SHV", 1, 6].item(), 0.1)
        self.assertEqual(multi.array.loc["2SHV", 6, 1].item(), 0.1)
        self.assertEqual(multi.array.loc["2SHV", 1, 3].item(), 0.2)

        # Creating another test case that is at analogous residues on the alignment but just shifted by 5
        multi.add(self.Dict_1ALI, "1ALI", seqList2)
        self.assertEqual(multi.array.loc["1ALI", 1, 6].item(), 0.1)
        self.assertEqual(multi.array.loc["1ALI", 6, 1].item(), 0.1)
        self.assertEqual(multi.array.loc["1ALI", 1, 3].item(), 0.2)

        multi.add(self.Dict_2SJR, "2SJR", seqList3)
        self.assertEqual(multi.array.loc["2SJR", 1, 6].item(), 0.1)
        self.assertEqual(multi.array.loc["2SJR", 6, 1].item(), 0.1)
        self.assertEqual(multi.array.loc["2SJR", 1, 3].item(), 0.2)

    def test_add_withNorm (self):

        args = Namespace(alignmentFile='tests/data/multi_net_test/PTP-KDY.fa', no_norm_struct=False)
        multi = MultiNetwork(args)
        structList = ["2SHV","1ALI","2SJR"]

        multi.array = xr.DataArray(
            0.0, 
            coords=dict(network=structList, firstResi=range(multi.size), secondResi=range(multi.size)), 
            dims=("network", "firstResi", "secondResi")
        )

        seqList1 = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
        seqList2 = [0, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]
        seqList3 = [0, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222]

        # Adds dict of dicts to MultiNetwork object and then tests each connection's weight
        multi.add(self.Dict_2SHV, "2SHV", seqList1)
        multi.add(self.Dict_1ALI, "1ALI", seqList2)
        multi.add(self.Dict_2SJR, "2SJR", seqList3)
        multi.normalizeStruct()

        self.assertEqual(multi.array.loc["2SHV", 1, 6].item(), 5)
        self.assertEqual(multi.array.loc["2SHV", 6, 1].item(), 5)
        self.assertEqual(multi.array.loc["2SHV", 1, 3].item(), 10)

        self.assertEqual(multi.array.loc["1ALI", 1, 6].item(), 5)
        self.assertEqual(multi.array.loc["1ALI", 6, 1].item(), 5)
        self.assertEqual(multi.array.loc["1ALI", 1, 3].item(), 10)

        self.assertEqual(multi.array.loc["2SJR", 1, 6].item(), 2.5)
        self.assertEqual(multi.array.loc["2SJR", 6, 1].item(), 2.5)
        self.assertEqual(multi.array.loc["2SJR", 1, 3].item(), 5)

    def test_NormalizeStruct (self):

        args = Namespace(alignmentFile='tests/data/multi_net_test/PTP-KDY.fa', no_norm_struct=False)
        multi = MultiNetwork(args)
        structList = ["2SHV","1ALI"]

        # Creates test array object
        multi.array = xr.DataArray(
            np.arange(18).reshape(2, 3, 3), 
            coords=dict(network=structList, firstResi=range(3), secondResi=range(3)), 
            dims=("network", "firstResi", "secondResi")
        )

        # Normalizes these values
        multi.normalizeStruct()

        self.assertEqual(multi.array.loc["2SHV", 1, 2].item(), (5/8) * 10)
        self.assertEqual(multi.array.loc["1ALI", 1, 2].item(), (14/17) * 10)

    def test_scaleMultiNet (self):

        args = Namespace(alignmentFile='tests/data/multi_net_test/PTP-KDY.fa', no_norm_struct=False, scale_value=20)
        multi = MultiNetwork(args)
        structList = ["2SHV","1ALI"]

        # Creates test array object
        multi.array = xr.DataArray(
            np.arange(18).reshape(2, 3, 3), 
            coords=dict(network=structList, firstResi=range(3), secondResi=range(3)), 
            dims=("network", "firstResi", "secondResi")
        )

        # Scales these values
        multi.scaleMultiNet()

        self.assertEqual(multi.array.loc["2SHV", 1, 2].item(), (5/17) * 20)
        self.assertEqual(multi.array.loc["1ALI", 1, 2].item(), (14/17) * 20)

    def test_sum (self):

        args = Namespace(alignmentFile='tests/data/multi_net_test/PTP-KDY.fa', no_norm_struct=True)
        multi = MultiNetwork(args)
        structList = ["2SHV","1ALI","2SJR"]

        multi.array = xr.DataArray(
            0.0, 
            coords=dict(network=structList, firstResi=range(multi.size), secondResi=range(multi.size)), 
            dims=("network", "firstResi", "secondResi")
        )

        seqList1 = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
        seqList2 = [0, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]
        seqList3 = [0, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222]

        # Adds dict of dicts to MultiNetwork object
        multi.add(self.Dict_2SHV, "2SHV", seqList1)
        multi.add(self.Dict_1ALI, "1ALI", seqList2)
        multi.add(self.Dict_2SJR, "2SJR", seqList3)

        # Calculates the sum across these three networks
        sumMulti = multi.sum()

        # Tests that positions at analogous positions are summed up
        self.assertEqual(round(sumMulti[1][3].item(), 1), 0.6)
        self.assertEqual(round(sumMulti[3][1].item(), 1), 0.6)
        self.assertEqual(round(sumMulti[1][6].item(), 1), 0.3)
        self.assertEqual(round(sumMulti[12][13].item(), 1), 0.4)

    def test_removeWeakEdges (self):

        args = Namespace(alignmentFile='tests/data/multi_net_test/PTP-KDY.fa', no_norm_struct=False, remove_weak_edges=20)
        multi = MultiNetwork(args)
        structList = ["2SHV","1ALI"]

        # Creates test array object
        # Two 5x5 arrays that have values from 0-50
        multi.array = xr.DataArray(
            np.arange(1,51,1).reshape(2, 5, 5), 
            coords=dict(network=structList, firstResi=range(5), secondResi=range(5)), 
            dims=("network", "firstResi", "secondResi")
        )

        # Gets original array's maximum value
        maxValue = multi.array.max().item()

        # Removes weak edges below 20% cutoff
        removedArray = multi.removeWeakEdges(multi.array)

        # Gets array of all values that are larger than zero but less than the cutoff
        # Ideally this should be empty
        cutoffArray = removedArray.where((removedArray < 0.2 * maxValue) & (removedArray > 0), drop=True)

        # Asserts that this array is empty (since there should be no values after removal that fit criteria)
        self.assertEqual(len(cutoffArray), 0)
 
if __name__ == '__main__':
    unittest.main()