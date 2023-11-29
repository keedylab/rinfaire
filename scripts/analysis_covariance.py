from multirin.analysis.Covariance import Covariance
from multirin.generate.MainFunctions import checkExtension
import argparse
import logging

def setupArguments ():

    # Creates argument parser object
    parser = argparse.ArgumentParser(
        description='Covariance',
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
        '-r',
        '--residue_pair', 
        default=False,
        action='store_true', 
        help="Calculates the covariance matrix based off of the residue-residue pairs across all networks"
    )

    parser.add_argument( 
        '--remove_weak_edges',
        type=int,
        help='Option to remove weak edges by a percent cutoff factor specified'
    )

    args = parser.parse_args()
    checkExtension(args.filename, '.pkl', "Input file must be in .pkl format")

    return args

def main ():
    
    # Sets up arguments inputted from user
    args = setupArguments()

    # Creates SumNetwork object and imports user inputted pickle file 
    covarianceObject = Covariance(args)
    covarianceObject.readPickle()

    # Calculates the covariance by residue pairs
    if args.residue_pair == True:
        covarianceObject.calculateCovarianceByResiPair()
        covarianceObject.calculateCorrelationByResiPair()

    # If the user wants to output a pickle file, then it calls on the outputPickle() function
    # If not then it gets visualized through PyVis which creates a .html output
    # Both files use the same file extension provided by the outputname argument
    if args.output_pickle == True:
        covarianceObject.exportPickle()
    else:
        covarianceObject.visualize()

if __name__ == "__main__":
    main()