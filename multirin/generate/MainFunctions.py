import argparse
import os
from .Structure import Structure
from .IndividualNetwork import IndividualNetwork
from .MultiNetwork import MultiNetwork

def setupArguments (multiFlag):

    # Creates argument parser object
    parser = argparse.ArgumentParser(
        description='Multiconformer model network generator',
        formatter_class=argparse.RawTextHelpFormatter
    )

    # Options that are specific to multi-model inputs
    if multiFlag == True:
    
        parser.add_argument(
            'structureFile', 
            help='List of input structures file in (.txt) format'
        )

        parser.add_argument(
            'alignmentFile', 
            help='Sequence alignment file in fasta (.fa) format'
        )

        parser.add_argument(
            '-m', 
            '--metadata',
            help='Table of metadata information to classify structures in (.csv) format'
        )

        parser.add_argument(
            '-f', 
            '--fetch', 
            action='store_true', 
            help='Fetch metadata from PDB API'
        )

        parser.add_argument( 
            '--no_norm_struct', 
            default=False,
            action='store_true', 
            help="Turns off the normalization of each structure's network in relation to the others"
        )

        parser.add_argument( 
            '--scale_multinet', 
            default=False,
            action='store_true', 
            help="Scales the MultiNetwork object to a range between 0 and 10"
        )

        parser.add_argument( 
            '--multinet_scale',
            default=10,
            type=int,
            help='Maximum edge weight the MultiNetwork should be scaled to'
        )

        parser.add_argument( 
            '--no_scale_sum_network', 
            default=False,
            action='store_true', 
            help="Turns off the scaling of the Sum Network"
        )

        parser.add_argument( 
            '--sum_network_scale',
            default=20,
            type=int,
            help='Maximum edge weight Sum Network should be scaled to'
        )

        parser.add_argument( 
            '--remove_weak_edges',
            type=int,
            help='Option to remove weak edges by a percent cutoff factor specified'
        )

        parser.add_argument( 
            '--output_info', 
            default=False,
            action='store_true', 
            help="Outputs summary statistics of the Multi Network"
        )
    
    # Options specific to single model inputs
    else:

        parser.add_argument(
            'structureFile', 
            help='Input structure file in (.pdb) format'
        )
    
    # Options that are shared between both programs

    parser.add_argument(
            'output',
            help='Output directory for all files generated'
        )

    parser.add_argument( 
        '--no_norm_resi', 
        default=False,
        action='store_true', 
        help="Turns off the normalization of each residue-residue pair in the network by the size of the residue"
    )

    parser.add_argument(
        '-a', 
        '--add_adjacent_residues', 
        default=False,
        action='store_true', 
        help="Adds residues spatially adjacent to the determined network"
    )

    args = parser.parse_args()

    # Checks of the file extension to ensure that it is the correct type
    # For multi-model inputs
    if multiFlag == True:
 
        checkExtension(args.structureFile, '.txt', "Structure file must be in .txt format")
        checkExtension(args.alignmentFile, '.fa', "Sequence file must be in .fa format")

        # Checks only if the user wants to use inputted labels
        if args.metadata != None:
            checkExtension(args.metadata, '.csv', "Labels file must be in .csv format")

    # For single model inputs
    else:
        checkExtension(args.structureFile, '.pdb', "Structure file must be in .pdb format")

    return args

def checkExtension (file, extension, errorMessage):
    fileSplit = os.path.splitext(file)
    fileExt = fileSplit[1]

    if fileExt != extension:
        raise argparse.ArgumentTypeError(errorMessage)
    
def readFile (multiFlag, args):
    
    fileList = []

    # Trivial case where you have a single model input
    if multiFlag == False:
        fileList.append(args.structureFile)
    
    # Case when you have a multi-model input
    else:
        # Opens file of pathnames
        inputStructureListFile = open(args.structureFile, "r")

        # Loops through each line (individual path) and adds it to list of pathnames
        for line in inputStructureListFile:
            fileList.append(line.rstrip())

        inputStructureListFile.close()
    
    return fileList

def generateIndividualNetworks (fileList, args):

    networkList = []

    # Loops over every pathname in the structure pathname list
    for structName in fileList:
        
        # Creates a structure object from the structure in pathname
        struct = Structure(structName, args)

        # Creates an individual network object from the structure object and then populates the network
        net = IndividualNetwork(struct, args)
        net.populateNetwork()

        if args.add_adjacent_residues == True:
            net.addAdjacentResidues()

        # Appends this to the list of networks
        networkList.append(net)

        # Creates the pyvis visualization (as an .html output)
        net.visualize()

    return networkList

def generateMultiNetwork (networkList, args):
    
    # Initializes an empty multi-network object 
    multi = MultiNetwork(args=args)

    # Adds networks from the list of individual networks
    multi.addNetworks(networkList)

    # Exports MultiNetwork object as a pickle file for further analysis
    multi.exportPickle()

    # Gets info about edges in MultiNetwork
    if args.output_info == True:
        multi.getInfo()

    # # Calculates the sum matrix of all the structures in the object
    # # Does this across each residue pairing
    # sumMatrix = multi.sum()
    # multi.visualize(sumMatrix, 'sumNetwork')

    return multi
    
