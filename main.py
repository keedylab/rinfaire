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

    args = parser.parse_args()

    print(args.structureListFile, args.alignmentFile, args.multi)
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

def generateNetworks (fileList, args):

    if args.multi == True:
        # Initializes an empty multi-network object 
        multi = MultiNetwork(args.alignmentFile)

    # Loops over every pathname in the structure pathname list
    for structPathName in fileList:
        
        # Creates a structure object from the structure in pathname
        struct = Structure(structPathName)

        # Creates an individual network object from the structure object and then populates the network
        net = IndividualNetwork(struct)
        net.populateNetwork()
        net.visualize()

        if args.multi == True:
            # Adds the individual network object to the multi-network object
            multi.Add(net.convertToAdjacency(), struct.getName()[0:4], struct.getFirstResi(), False)
            print(multi.array)

    if args.multi == True:
        return multi
    else:
        return None

def main ():
    args = setupArguments()
    fileList = readFileList(args)
    multiNetwork = generateNetworks(fileList, args)

if __name__ == "__main__":
    main()