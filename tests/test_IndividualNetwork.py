from multirin import Structure, IndividualNetwork
from argparse import Namespace
import gemmi
import unittest

# Define class to test the program
class testIndividualNetwork (unittest.TestCase):

    @classmethod
    def setUpClass (self):
        self.struct1 = Structure.Structure('tests/data/synth_structs/1111_Test.pdb', None)
    
    # Function to test alt-conf finder function
    def test_findAltConfs (self):

        args1 = Namespace(no_norm_resi=False)
        net1 = IndividualNetwork(self.struct1, args1)

        # Tests whether the correct residues are being flagged as having alt-confs
        self.assertEqual(list(net1.findAltConfAtoms().keys()), [1])
        # Tests whether the number of alt-conf atoms is correct
        self.assertEqual(len(net1.findAltConfAtoms()[1]), 28)

    # TODO:
    # def test_flagAmideHydrogenOnlyResidues (self):

    def test_updateEdge_withNorm (self):

        args1 = Namespace(no_norm_resi=False)
        net1 = IndividualNetwork(self.struct1, args1)

        # Sets residue-residue counter
        # Then updates edges (twice for 1,3 and once for 1,5)
        counter = 0
        counter = net1.updateEdge(1, 3, counter)
        counter = net1.updateEdge(1, 3, counter)
        counter = net1.updateEdge(1, 5, counter)
        # Ensures that there are only two resi-resi connections (since one edge had two atoms added and therefore "duplicated")
        self.assertEqual(counter, 2)

        # Checks the weight of the connections after residue-level normalization
        # First connection has two atom connections – therefore weight is 2 (times 10) divided by number of atoms
        # Second connection only has one atom connection
        self.assertEqual(net1.network[1][3]['weight'], (20 / (28 + 24)))
        self.assertEqual(net1.network[1][5]['weight'], (10 / (28 + 19)))

    def test_updateEdge_withoutNorm (self):

        args1 = Namespace(no_norm_resi=True)
        net1 = IndividualNetwork(self.struct1, args1)

        # Sets residue-residue counter
        # Then updates edges (twice for 1,3 and once for 1,5)
        counter = 0
        counter = net1.updateEdge(1, 3, counter)
        counter = net1.updateEdge(1, 3, counter)
        counter = net1.updateEdge(1, 5, counter)
        # Ensures that there are only two resi-resi connections (since one edge had two atoms added and therefore "duplicated")
        self.assertEqual(counter, 2)

        # Same check but without residue-level normalization
        self.assertEqual(net1.network[1][3]['weight'], 0.2)
        self.assertEqual(net1.network[1][5]['weight'], 0.1)

    # TODO:
    # def test_populateNetwork (self):

    def test_convertToAdjacency (self):

        # Creates network that's the same as previous functions
        args1 = Namespace(no_norm_resi=True)
        net1 = IndividualNetwork(self.struct1, args1)
        counter = 0
        counter = net1.updateEdge(1, 3, counter)
        counter = net1.updateEdge(1, 3, counter)
        counter = net1.updateEdge(1, 5, counter)

        # Creates the correct dict representation to compare the output with
        compareDict = {
            1: {3: {'weight': 0.2}, 5: {'weight': 0.1}},
            3: {1: {'weight': 0.2}},
            5: {1: {'weight': 0.1}}
        }

        self.assertEqual(net1.convertToAdjacency(), compareDict)

if __name__ == '__main__':
    unittest.main()