from multirin.analysis.ResiduesOfInterest import ResiduesOfInterest
from multirin.generate.MainFunctions import checkExtension
import argparse
import logging

def setupArguments ():

    # Creates argument parser object
    parser = argparse.ArgumentParser(
        description='Residues of Interest',
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(
        'filename', 
        help='Pickle (.pkl) file of the SumNetwork object' 
    )

    parser.add_argument(
        'outputname', 
        help='Name of the output file' 
    )

    parser.add_argument(
        '-i',
        '--input_set', 
        help='Name of the file with the set of residues to input and compare' 
    )

    parser.add_argument(
        '-c',
        '--col',
        default=0,
        help='Column of input csv to use' 
    )

    parser.add_argument(
        '-a',
        '--include_adjacent_residues', 
        help='Flag to check adjacent residues to network to check for hits. Must include reference structure file.' 
    )

    parser.add_argument(
        '-s',
        '--find_significance', 
        help='Flag to find the significance of the overlap of the input set compared to all residues. Must include reference structure file.' 
    )

    parser.add_argument(
        '--n_iter_sig_test',
        default=100,
        help='Number of iterations to perform the significance test' 
    )

    parser.add_argument(
        '--no_normalize_by_total', 
        default=False,
        action='store_true',
        help='To be used with find_significance. Does not normalize number of adjacent network residues by total number of residues.' 
    )

    parser.add_argument(
        '--histogram', 
        default=False,
        action='store_true',
        help='To be used with find_significance. Outputs traditional histogram of the distributions.' 
    )
    
    parser.add_argument(
        '--cumulative_histogram', 
        default=False,
        action='store_true',
        help='To be used with find_significance. Outputs the cumulative instead of traditional histogram.' 
    )

    parser.add_argument(
        '-p',
        '--output_pickle',
        default=False,
        action='store_true',
        help='Flag to output pickle file instead of .html visualization' 
    )

    args = parser.parse_args()
    checkExtension(args.filename, '.pkl', "Input file must be in .pkl format")

    return args

def main ():
    
    # Sets up arguments inputted from user
    args = setupArguments()

    # Creates ResiduesOfInterest object and imports user inputted pickle file 
    resiObject = ResiduesOfInterest(args)
    resiObject.readPickle()

    # Reads in the input set of residues from file
    resiObject.readInputSetFile()
    resiObject.findOverlapInputSet()
    
    if args.find_significance != None:
        resiObject.findSignificance()

    resiObject.labelGraphOverlap()
    resiObject.visualize(resiObject.overlapGraph, "overlapGraph")
    resiObject.visualize(resiObject.adjResisNetwork.network, "adjacentResiNetwork")


    # # If the user wants to output a pickle file, then it calls on the outputPickle() function
    # # If not then it gets visualized through PyVis which creates a .html output
    # # Both files use the same file extension provided by the outputname argument
    # if args.output_pickle == True:
    #     sumNetObject.exportPickle()
    # else:
    #     sumNetObject.visualize()

if __name__ == "__main__":
    main()