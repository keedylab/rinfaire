from multirin.generate.Structure import Structure 
from argparse import Namespace
import gemmi
import unittest

# Define class to test the program
class testStructure (unittest.TestCase):

    @classmethod
    def setUpClass (self):

        # structureListFile='tests/data/synth_structs/TestInputFiles.txt'
        self.args = Namespace()
        self.struct1 = Structure('tests/data/synth_structs/1111_Test.pdb', self.args)
    
    # Function to test model of structure
    def test_model (self):
        self.assertEqual(type(self.struct1.model), gemmi.Structure)

    # Function to test the name of structure
    def test_name (self):
        self.assertEqual(self.struct1.name, '1111_Test')

    # Function to test the sequence 
    def test_sequence (self):
        self.assertEqual(len(self.struct1.sequence), 9)
        self.assertEqual(self.struct1.sequence[5]['name'], 'LEU')
        self.assertEqual(self.struct1.sequence[5]['atomcount'], 19)

    def test_sequencelist (self):
        self.assertEqual(len(self.struct1.sequenceList), 10)

    def test_firstresi (self):
        self.assertEqual(self.struct1.getFirstResi(), 1)
 
if __name__ == '__main__':
    unittest.main()