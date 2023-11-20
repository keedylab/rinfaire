from multirin.analysis.SumNetwork import SumNetwork
from multirin.generate.MainFunctions import checkExtension
import argparse
import logging

def setupArguments ():

    # Creates argument parser object
    parser = argparse.ArgumentParser(
        description='SumNetwork',
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(
        'filename', 
        help='Pickle (.pkl) file of the MultiNetwork object' 
    )

    parser.add_argument(
        'outputname', 
        help='Name of the output file' 
    )

    parser.add_argument(
        '-p',
        '--output_pickle',
        default=False,
        action='store_true',
        help='Flag to output pickle file instead of .html visualization' 
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
        '--no_resize_by_degree', 
        default=False,
        action='store_true', 
        help="Turns off node resizing by degree"
    )

    parser.add_argument( 
        '--resize_by_degree_scale',
        default=1,
        type=int,
        help='Maximum edge weight Sum Network should be scaled to'
    )

    parser.add_argument( 
        '--remove_subgraphs',
        default=5,
        type=int,
        help='Removes subgraphs that have less than specified number of residues. To turn off removal set this to 0'
    )

    parser.add_argument(
        '-r',
        '--seq_to_ref',
        help='Sets a reference sequence for all the numbering'
    )

    parser.add_argument( 
        '-c',
        '--detect_communities', 
        default=False,
        action='store_true', 
        help="Detects communities within the sum graph"
    )

    args = parser.parse_args()
    checkExtension(args.filename, '.pkl', "Input file must be in .pkl format")

    return args

def main ():
    
    # Sets up arguments inputted from user
    args = setupArguments()

    # Creates SumNetwork object and imports user inputted pickle file 
    sumNetObject = SumNetwork(args)
    sumNetObject.readPickle()

    # Calculates the sum network and performs scaling and removal of weak edges if flags are specified
    sumNetObject.calculateSum()

    # If the user wants to output a pickle file, then it calls on the outputPickle() function
    # If not then it gets visualized through PyVis which creates a .html output
    # Both files use the same file extension provided by the outputname argument
    if args.output_pickle == True:
        sumNetObject.exportPickle()
    else:
        sumNetObject.visualize()

if __name__ == "__main__":
    main()