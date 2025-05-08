from multirin.analysis.Covariance import Covariance
from multirin.generate.MainFunctions import checkExtension
import argparse
import logging

def setupArguments ():

    # Creates argument parser object
    parser = argparse.ArgumentParser(
        description='PCA',
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

    covarianceObject.calculateCovarianceByResiPair(scaleFlag=True)
    covarianceObject.runPCA()


    

if __name__ == "__main__":
    main()