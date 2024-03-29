from multirin.generate.MultiNetwork import MultiNetwork
from multirin.generate.Subset import readPickle, exportPickle, generateSubsets
from multirin.generate.MainFunctions import checkExtension
import argparse
import logging

def setupArguments ():

    # Creates argument parser object
    parser = argparse.ArgumentParser(
        description='Subset',
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(
        'filename', 
        help='Pickle (.pkl) file of the MultiNetwork object' 
    )

    parser.add_argument(
        'subset', 
        help='Name of the metadata column to subset by' 
    )

    parser.add_argument(
        'outputname', 
        help='Name of the output file' 
    )

    parser.add_argument(
        '-g',
        '--group',
        help='Specific group within the column to select' 
    )

    parser.add_argument(
        '-d',
        '--make_discrete',
        type=float,
        default=None,
        nargs=3,
        help='Makes continuous numerical values in column into discrete bins. Must specify as arguments (in order): bin size, starting bin, ending bin. \n Ex. for resolutions between 1 and 3A with bins of 0.5A: -d 0.5 1.0 3.0' 
    )

    args = parser.parse_args()
    checkExtension(args.filename, '.pkl', "Input file must be in .pkl format")

    return args

def main ():
    
    # Sets up arguments inputted from user
    args = setupArguments()

    # Creates MultiNetwork object
    multinet = readPickle(args.filename)

    # Generates subset MultiNetworks by subsetting by the classifier
    if args.group != None:
        subsetMultiNetworks = generateSubsets(multinet, args.subset, groupName=args.group, makeDiscreteValue=args.make_discrete)
    else:
        subsetMultiNetworks = generateSubsets(multinet, args.subset, makeDiscreteValue=args.make_discrete)

    # Exports a series of pickle files
    exportPickle(subsetMultiNetworks, args.outputname)

if __name__ == "__main__":
    main()