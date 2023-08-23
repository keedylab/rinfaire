from multirin import MultiNetwork 
from argparse import Namespace
import unittest
import numpy as np

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
            1: {3: {'weight': 0.2}, 5: {'weight': 0.1}},
            3: {1: {'weight': 0.2}},
            5: {1: {'weight': 0.1}},
            10: {11: {'weight': 0.4}},
            11: {10: {'weight': 0.4}}
        }

    # Function to test residue conversion function
    def test_oneToAll (self):

        args = Namespace(alignmentFile='tests/data/multi_net_test/PTP-KDY.fa', no_normalizing_struct=False)
        multi = MultiNetwork(args)

        # Tests with inputs of:
        #   Fake PDB entry: 2SHV, 3LJR
        #   Starting residue: 5, 200
        #   Residue to convert: 10, 210
        result1 = multi.oneToAll('2SHV', 5, 10)
        result2 = multi.oneToAll('3LJR', 200, 210)
        self.assertEqual(result1, 7)
        self.assertEqual(result2, 14)

    def test_add_withoutNorm (self):

        args = Namespace(alignmentFile='tests/data/multi_net_test/PTP-KDY.fa', no_normalizing_struct=True)
        multi = MultiNetwork(args)

        multi.add(self.Dict_2SHV, "2SHV", 1)
        self.assertEqual(multi.array[0][1][6], 0.1)
        self.assertEqual(multi.array[0][6][1], 0.1)
        self.assertEqual(multi.array[0][1][3], 0.2)

        # Creating another test case that is at analogous residues on the alignment but just shifted differently
        multi.add(self.Dict_1ALI, "1ALI", 5)
        self.assertEqual(multi.array[1][1][6], 0.1)
        self.assertEqual(multi.array[1][6][1], 0.1)
        self.assertEqual(multi.array[1][1][3], 0.2)

    def test_add_withNorm (self):

        args = Namespace(alignmentFile='tests/data/multi_net_test/PTP-KDY.fa', no_normalizing_struct=False)
        multi = MultiNetwork(args)

        # Adds dict of dicts to MultiNetwork object and then tests each connection's weight
        multi.add(self.Dict_2SHV, "2SHV", 1)
        self.assertEqual(multi.array[0][1][6], 5)
        self.assertEqual(multi.array[0][6][1], 5)
        self.assertEqual(multi.array[0][1][3], 10)

        multi.add(self.Dict_1ALI, "1ALI", 5)
        self.assertEqual(multi.array[1][1][6], 5)
        self.assertEqual(multi.array[1][6][1], 5)
        self.assertEqual(multi.array[1][1][3], 10)

    def test_NormalizeStruct (self):

        args = Namespace(alignmentFile='tests/data/multi_net_test/PTP-KDY.fa', no_normalizing_struct=False)
        multi = MultiNetwork(args)

        # Creates a simple testing array with values
        testArray = np.array(
            [[0, 0, 0, 0],
             [0, 0, 0.5, 0.3],
             [0, 0.5, 0, 0.2],
             [0, 0.3, 0.2, 0]]
        )

        # Creates a normalized array from the testing array using the normalizeStruct function
        normArray = multi.normalizeStruct(testArray)
        self.assertEqual(normArray[1][2], 10)
        self.assertEqual(normArray[2][1], 10)
        self.assertEqual(normArray[1][3], 6)
        self.assertEqual(normArray[2][3], 4)

    def test_sum (self):

        args = Namespace(alignmentFile='tests/data/multi_net_test/PTP-KDY.fa', no_normalizing_struct=True)
        multi = MultiNetwork(args)

        # Adds dict of dicts to MultiNetwork object
        multi.add(self.Dict_2SHV, "2SHV", 1)
        multi.add(self.Dict_1ALI, "1ALI", 5)
        multi.add(self.Dict_2SJR, "2SJR", 1)

        # Calculates the sum across these three networks
        sumMulti = multi.sum()

        # Tests that positions at analogous positions are summed up
        self.assertEqual(round(sumMulti[1][3], 1), 0.6)
        self.assertEqual(round(sumMulti[3][1], 1), 0.6)
        self.assertEqual(round(sumMulti[1][6], 1), 0.3)
        self.assertEqual(round(sumMulti[12][13], 1), 0.4)
 
if __name__ == '__main__':
    unittest.main()