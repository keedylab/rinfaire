from MultiNetwork import MultiNetwork
import unittest

# Define class to test the program
class testGenerateMultiNetwork(unittest.TestCase):

    # Function to conversion function
    def test_oneToAll(self):

        multi = MultiNetwork()

        result = MultiNetwork.oneToAll('6B8Z', )
        #self.assertEqual(result, 4)
 
if __name__ == '__main__':
    unittest.main()