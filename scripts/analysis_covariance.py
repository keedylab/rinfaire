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
        'outputdir', 
        help='Name of the output directory' 
    )

    parser.add_argument( 
        '-c',
        '--correlation', 
        default=False,
        action='store_true', 
        help="Calculates Pearson correlation coefficient matrix instead of covariance matrix"
    )

    parser.add_argument(
        '-p',
        '--output_pickle',
        default=False,
        action='store_true',
        help='Flag to output pickle file instead of .html visualization' 
    )

    parser.add_argument( 
        '-s',
        '--structure_pair', 
        default=False,
        action='store_true', 
        help="Calculates the covariance matrix based off of the structure-structure pairs across all residues (instead of default where it calculates residue-residue pairs across all structures)"
    )

    parser.add_argument( 
        '--remove_weak_edges',
        type=int,
        help='Option to remove weak edges by a percent cutoff factor specified'
    )

    parser.add_argument( 
        '-l',
        '--cluster_corr', 
        default=False,
        action='store_true', 
        help="Clusters the correlation matrix using heirarchical clustering"
    )

    parser.add_argument(
        '--vis_from_sum', 
        help='Maps the covariance/correlation values onto the SumNetwork (must also specify .pkl file of SumNetwork object)' 
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
    if args.structure_pair == False:

        # Calculates the correlation matrix
        if args.correlation == True:
            covarianceObject.calculateCorrelationByResiPair()
        
        # Calculates the covariance matrix
        else:
            covarianceObject.calculateCovarianceByResiPair(scaleFlag=False)

    if args.cluster_corr == True:
        covarianceObject.clusterCorrMatrix()

    # If the user wants to output a pickle file, then it calls on the outputPickle() function
    # If the user wants to map either matrix to the SumNetwork then it calls on visualizeFromSumNetwork
    # If not then it gets visualized through PyVis which creates a .html output
    # Both files use the same file extension provided by the outputname argument
    if args.output_pickle == True:
        covarianceObject.exportPickle()

    elif args.vis_from_sum == True:
        if args.correlation == True:
            covarianceObject.visualizeFromSumNetwork('correlation')
        else:
            covarianceObject.visualizeFromSumNetwork('covariance')

    else:
        if args.correlation == True:
            covarianceObject.visualizeMatrix('correlation')
        else:
            covarianceObject.visualizeMatrix('covariance')

if __name__ == "__main__":
    main()