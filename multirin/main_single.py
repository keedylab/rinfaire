from IndividualNetwork import IndividualNetwork
from Structure import Structure, Chain
import argparse
import os

def setupArguments ():

    # Creates argument parser object
    parser = argparse.ArgumentParser(
        description='Multiconformer model network generator',
        formatter_class=argparse.RawTextHelpFormatter
    )
    
    # Adds individual arguments
    parser.add_argument(
        'structureFile', 
        help='Input structure file in (.pdb) format'
    )

    parser.add_argument( 
        '--no_norm_resi', 
        default=False,
        action='store_true', 
        help="Turns off the normalization of each residue-residue pair in the network by the size of the residue"
    )

    args = parser.parse_args()

    # Checks the extension of the file given to ensure that it is the correct type
    checkExtension(args.structureFile, '.pdb', "Structure file must be in .pdb format")
    
    return args

def checkExtension (file, extension, errorMessage):
    fileSplit = os.path.splitext(file)
    fileExt = fileSplit[1]

    if fileExt != extension:
        raise argparse.ArgumentTypeError(errorMessage)

def generateNetwork (args):

    # Creates a structure object from the structure in pathname, then gets the chains into a list
    struct = Structure(args.structureFile, args)
    chainList = struct.getChains()

    print(struct)
    print(chainList[0].structure)

    # Creates list of networks
    netList = []

    # TODO: Maybe create a dictionary for the chains, then iterate through dict to add all the chains in

    # # Iterates over each chain in the list and creates a network for each
    # for chain in chainList:

    #     # Creates an individual network object from the structure object and then populates the network
    #     net = IndividualNetwork(chain, args)
    #     net.populateNetwork()

    #     # Appends this to the list of networks
    #     netList.append(net)

    #     # Creates the pyvis visualization (as an .html output)
    #     net.visualize()

    return netList

def main ():
    args = setupArguments()
    networkList = generateNetwork(args)

if __name__ == "__main__":
    main()