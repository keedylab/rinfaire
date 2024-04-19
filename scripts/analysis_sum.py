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
        '-s',
        '--subset',
        nargs = 2,
        help='Specifies to calculate sum based on a specific subset of the metadata. Takes in two arguments, the first is the column to subset by and the second is the group of interest within that subset. (note: MUST have included metadata to MultiNetwork)'
    )

    parser.add_argument( 
        '--subset_all',
        help='Same as subset but does this for all groups in the column. (note: MUST have included metadata to MultiNetwork)'
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
        default=0,
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

    parser.add_argument( 
        '--n_communities',
        type=int,
        help='Number of communities to be detected (must be used with -c flag)'
    )

    parser.add_argument( 
        '--output_modularity', 
        default=False,
        action='store_true', 
        help="Outputs .csv file with modularity values for every value of k (must be used with -c flag)"
    )

    parser.add_argument( 
        '--output_graph_info', 
        default=False,
        action='store_true', 
        help="Outputs information about the graph including a .csv file with degrees of every node"
    )

    parser.add_argument( 
        '--mst', 
        default=False,
        action='store_true', 
        help="Outputs Maximum Spanning Tree instead of regular graph"
    )

    parser.add_argument( 
        '--keep_nan', 
        default=False,
        action='store_true', 
        help="Used with --seq_to_ref. Keeps nodes that don't shift to any residue in input structure."
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
    if args.subset != None:
        sumNetObject.generateSumNetworkSubset(args.subset[0], groupName=args.subset[1])
    elif args.subset_all != None:
        sumNetObject.generateSumNetworkSubset(args.subset_all)
    else:
        sumNetObject.generateSumNetworkAll()

    # Constructs Maximum Spanning Tree as graph output if specified, if not then just outputs normal graph
    if args.mst == True:
        sumNetObject.constructMaxSpanningTrees()

    else:
        # Creates the graph representation in networkX
        sumNetObject.constructGraphs()

    # If the user wants information about the graph to be printed
    if args.output_graph_info == True:
        sumNetObject.exportGraphInfo()

    # If the user wants to output a pickle file, then it calls on the outputPickle() function
    # If not then it gets visualized through PyVis which creates a .html output
    # Both files use the same file extension provided by the outputname argument
    if args.output_pickle == True:
        sumNetObject.exportPickle()
    else:
        sumNetObject.visualize()

if __name__ == "__main__":
    main()