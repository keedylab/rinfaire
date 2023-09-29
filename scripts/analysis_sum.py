from multirin.analysis.SumNetwork import SumNetwork
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

    args = parser.parse_args()

    return args

def main ():
    
    # Sets up arguments inputted from user
    args = setupArguments()

    # Creates SumNetwork object and imports user inputted pickle file 
    sumNetObject = SumNetwork(args)

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