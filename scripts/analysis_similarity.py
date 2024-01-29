from multirin.analysis.Similarity import Similarity
from multirin.generate.MainFunctions import checkExtension
import argparse
import logging

def setupArguments ():

    # Creates argument parser object
    parser = argparse.ArgumentParser(
        description='Similarity',
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(
        'filename', 
        help='Pickle (.pkl) file of the MultiNetwork object' 
    )

    parser.add_argument(
        'outputdir', 
        help='Name of the output directory' 
    )

    parser.add_argument(
        '-p',
        '--output_pickle',
        default=False,
        action='store_true',
        help='Flag to output pickle file instead of .html visualization' 
    )

    parser.add_argument(
        '-d',
        '--distance_metric',
        default='jaccard',
        help='Flag to change the distance function between graphs. Options are Jaccard, DeltaCon, ...' 
    )

    args = parser.parse_args()
    checkExtension(args.filename, '.pkl', "Input file must be in .pkl format")

    return args

def main ():
    
    # Sets up arguments inputted from user
    args = setupArguments()

    # Creates Similarity class object and imports user inputted pickle file 
    similarityObject = Similarity(args)
    similarityObject.readPickle()

    similarityObject.distanceMatrix()
    similarityObject.heirClustering()
    similarityObject.visualizeMatrix()
    

if __name__ == "__main__":
    main()