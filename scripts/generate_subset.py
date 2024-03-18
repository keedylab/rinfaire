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

    args = parser.parse_args()
    checkExtension(args.filename, '.pkl', "Input file must be in .pkl format")

    return args

def main ():
    
    # Sets up arguments inputted from user
    args = setupArguments()

    # Creates MultiNetwork object
    multinet = readPickle(args.filename)

    # Generates subset arrays by subsetting by the classifier
    if args.group != None:
        subsetArrays = generateSubsets(multinet, args.subset, groupName=args.group)
    else:
        subsetArrays = generateSubsets(multinet, args.subset)

    # Exports a series of pickle files
    exportPickle(subsetArrays, args.outputname)

if __name__ == "__main__":
    main()