from multirin import main, MultiNetwork 
from argparse import Namespace
import unittest

# Define class to test the program
class testGenerateMultiNetwork(unittest.TestCase):

    @classmethod
    def setUpClass(self):

        self.args = Namespace(alignmentFile='tests/data/multi_net_test/PTP-KDY.fa')
        self.multi = MultiNetwork(self.args)

    # Function to test residue conversion function
    def test_oneToAll(self):

        # Tests with inputs of:
        #   Fake PDB entry: 2SHV, 3LJR
        #   Starting residue: 5, 200
        #   Residue to convert: 10, 210
        result1 = self.multi.oneToAll('2SHV', 5, 10)
        result2 = self.multi.oneToAll('3LJR', 200, 210)
        self.assertEqual(result1, 7)
        self.assertEqual(result2, 14)
 
if __name__ == '__main__':
    unittest.main()