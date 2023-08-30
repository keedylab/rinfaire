from IndividualNetwork import IndividualNetwork
from Structure import Structure
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
      
    # Creates a structure object from the structure in pathname
    struct = Structure(args.structureFile, args)

    # Creates an individual network object from the structure object and then populates the network
    net = IndividualNetwork(struct, args)
    net.populateNetwork()
    net.visualize()

    return net

def main ():
    args = setupArguments()
    network = generateNetwork(args)

if __name__ == "__main__":
    main()