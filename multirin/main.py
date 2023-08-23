from IndividualNetwork import IndividualNetwork
from Structure import Structure
from MultiNetwork import MultiNetwork
import argparse

def setupArguments ():

    # Creates argument parser object
    parser = argparse.ArgumentParser(description='Multiconformer model network generator')
    
    # Adds individual arguments
    parser.add_argument(
        'structureListFile', 
        help='Text file (.txt) of the pathnames to all of the structures'
    )

    parser.add_argument(
        'alignmentFile', 
        help='Sequence alignment file in fasta (.fa) format'
    )

    parser.add_argument(
        '-m', 
        '--multi', 
        action='store_true', 
        help='Flag that if present combines and analyzes the multiple networks inputted together'
    )

    parser.add_argument( 
        '--no_normalizing_resi', 
        default=False,
        action='store_true', 
        help="Turns off the normalization of each residue-residue pair in the network by the size of the residue"
    )

    parser.add_argument( 
        '--no_normalizing_struct', 
        default=False,
        action='store_true', 
        help="Turns off the normalization of each structure's network in relation to the others"
    )

    args = parser.parse_args()
    return args

def readFileList (args):
    
    fileList = []

    # Opens file of pathnames
    inputStructureListFile = open(args.structureListFile, "r")

    # Loops through each line (individual path) and adds it to list of pathnames
    for line in inputStructureListFile:
        fileList.append(line.rstrip())

    inputStructureListFile.close()
    return fileList

def generateIndividualNetworks (fileList, args):

    networkList = []

    # Loops over every pathname in the structure pathname list
    for structPathName in fileList:
        
        # Creates a structure object from the structure in pathname
        struct = Structure(structPathName, args)

        # Creates an individual network object from the structure object and then populates the network
        net = IndividualNetwork(struct, args)
        net.populateNetwork()
        #net.visualize()

        # Appends the network generated into the list of networks
        networkList.append(net)

    return networkList
    
def generateMultiNetwork(networkList, args):
    
    # Initializes an empty multi-network object 
    multi = MultiNetwork(args)

    # Loops over each network in the list
    for net in networkList:

        # Adds the individual network object to the multi-network object
        multi.add(
            net.convertToAdjacency(), 
            net.struct.getName()[0:4], 
            net.struct.getFirstResi()
        )

    # Calculates the sum matrix of all the structures in the object
    # Does this across each residue pairing
    sumMatrix = multi.sum()
    multi.visualize(sumMatrix, 'sumNetwork')

    return multi

def main ():
    args = setupArguments()
    fileList = readFileList(args)
    networkList = generateIndividualNetworks(fileList, args)
    
    if args.multi == True:
        multi = generateMultiNetwork(networkList, args)

if __name__ == "__main__":
    main()