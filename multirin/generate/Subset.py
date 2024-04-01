import xarray as xr
import numpy as np
import pandas as pd
import pickle
from multirin.generate.MultiNetwork import MultiNetwork

def readPickle (inputFile):
    
    """
    Function that opens MultiNetwork pickle file and returns it as object
    """
        
    # Opens pickle file
    with open(inputFile, 'rb') as pickleFile:
        multinet = pickle.load(pickleFile)

    if multinet.metadata is None:
        raise FileNotFoundError('Must have metadata associated with MultiNetwork!')

    return multinet

def generateSubsets (multinet, classifier, groupName=None, makeDiscreteValue=None):

    """
    Function that takes in a MultiNetwork object and splits them into subsets based on the classifier (column from metadata)
    
    Inputs: 
    multinet – MultiNetwork object
    classifier – Label of column in the MultiNetwork metadata .csv that you want to group by
    groupName – Specific group within the column that is of interest (optional argument)
    makeDiscreteValue – List of parameters specifying how you want to split up a continuous variable into discrete bins (optional argument)

    Output:
    subsetArrays – Dictionary of subset arrays with keys as groups
    """

    # Initializes dictionary of arrays that hold each subset array
    subsetMultiNetworks = {}

    # Condition when there is a specified group within the classifier column that is of interest
    if groupName != None:

        # Condition where we need to split a continuous variable in metadata column into discrete bins
        if makeDiscreteValue != None:
            
            # Groups by column after first using pandas' cut function to split column into discrete bins
            # makeDiscreteValue list is the user specified parameters for bins supplied by the make_discrete flag. Element 0 is interval, 1 is lower bound, and 2 is upper bound
            groupNameInterval = pd.Interval(left=groupName.split('-')[0], right=groupName.split('-')[1], closed='right')
            dfGrouped = multinet.metadata.groupby(pd.cut(multinet.metadata[classifier], np.arange(makeDiscreteValue[1], makeDiscreteValue[2] + makeDiscreteValue[0], makeDiscreteValue[0])), observed=True).get_group(groupNameInterval)
        
        else:
            # Groups by the subset classifier and then creates groups, then gets specific group
            dfGrouped = multinet.metadata.groupby(by=classifier, observed=True).get_group(groupName)
        
        # Creates dictionary of group's name and dataframe associated with it 
        groups = {groupName: dfGrouped}

    # Otherwise this will generate subsets for all instances in that classifier column
    else:

        # Condition where we need to split a continuous variable in metadata column into discrete bins
        if makeDiscreteValue != None:

            # Groups by column after first using pandas' cut function to split column into discrete bins
            # makeDiscreteValue list is the user specified parameters for bins supplied by the make_discrete flag. Element 0 is interval, 1 is lower bound, and 2 is upper bound
            dfGrouped = multinet.metadata.groupby(pd.cut(multinet.metadata[classifier], np.arange(makeDiscreteValue[1], makeDiscreteValue[2] + makeDiscreteValue[0], makeDiscreteValue[0])), observed=True)

        else:
            dfGrouped = multinet.metadata.groupby(by=classifier, observed=True)

        groups = dict(list(dfGrouped))

    # Gets list of structures in the MultiNetwork
    structsInMultiNet = set(multinet.array.get_index('network').to_list())

    # Iterates over each group
    for group in groups:

        # Gets list of structures in the group, then finds intersection of this with MultiNetwork structures
        structsInGroup = groups[group]['ID'].to_list()
        interStructs = list(set(structsInMultiNet).intersection(set(structsInGroup)))

        # Checks if intersection is empty or not
        if interStructs != []:
            
            # Selects subset of structures in group from the MultiNetwork object to create a smaller MultiNetwork array
            subsetArray = multinet.array.loc[interStructs, :, :]

            # Selects subset of metadata
            subsetMetadata = multinet.metadata.loc[multinet.metadata['ID'].isin(interStructs)]

            # Selects subset of sequence alignment
            subsetSeqAln = {k: multinet.seqaln[k] for k in interStructs}

            # Combines the subset array, metadata, and sequence alignment into a new MultiNetwork object
            # Puts this in a dictionary of MultiNetworks
            subsetMultiNetworks[group] = MultiNetwork(array=subsetArray, seqaln=subsetSeqAln, metadata=subsetMetadata)

        else:
            print(f"No structures in MultiNetwork for group: {group}")

    return subsetMultiNetworks

def exportPickle (subsetMultiNetworks, outputname):

    """
    Function that creates a pickle file for each subset array
    """

    # Loops through each subset array
    for subset in subsetMultiNetworks:

        # String formatting for proper output file name
        subsetString = str(subset).replace(', ','-').replace('(','').replace(']','')

        # Creates new pickle (.pkl) file and then dumps the entire class object into the pickle file
        with open(f'{outputname}_{subsetString}.pkl', 'wb') as pickleFile:
            pickle.dump(subsetMultiNetworks[subset], pickleFile)