import xarray as xr
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

    if isinstance(multinet.metadata, pd.DataFrame) == False:
        print("Error: Must have metadata associated with MultiNetwork!")
        return

    return multinet

def generateSubsets (multinet, classifier, groupName=None):

    """
    Function that takes in a MultiNetwork object and splits them into subsets based on the classifier (column from metadata)
    
    Inputs: 
    multinet – MultiNetwork object
    classifier – Label of column in the MultiNetwork metadata .csv that you want to group by
    groupName – Specific group within the column that is of interest (optional argument)

    Output:
    subsetArrays – Dictionary of subset arrays with keys as groups
    """

    # Initializes dictionary of arrays that hold each subset array
    subsetMultiNetworks = {}

    # Condition when there is a specified group within the classifier column that is of interest
    if groupName != None:

        # Groups by the subset classifier and then creates groups, then gets specific group
        dfGrouped = multinet.metadata.groupby(by=classifier).get_group(groupName)
        
        # Creates dictionary of group's name and dataframe associated with it 
        groups = {groupName: dfGrouped}

    # Otherwise this will generate subsets for all instances in that classifier column
    else:
        dfGrouped = multinet.metadata.groupby(by=classifier)
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

def exportPickle (subsetArrays, outputname):

    """
    Function that creates a pickle file for each subset array
    """

    # Loops through each subset array
    for subset in subsetArrays:

        # Creates new pickle (.pkl) file and then dumps the entire class object into the pickle file
        with open(f'{outputname}_{subset}.pkl', 'wb') as pickleFile:
            pickle.dump(subsetArrays[subset], pickleFile)