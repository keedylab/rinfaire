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
    
    args = setupArguments()
    sumNetObject = SumNetwork(args)

if __name__ == "__main__":
    main()