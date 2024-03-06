import xarray as xr
import pandas as pd

def generateSubsets (multinet, classifier):

    """
    Function that takes in a MultiNetwork object and splits them into subsets based on the classifier (column from metadata)
    
    Inputs: 
    multinet – MultiNetwork object
    classifier – Label of column in the MultiNetwork metadata .csv that you want to group by

    Output:
    subsetArrays – Dictionary of subset arrays with keys as groups
    """

    # Initializes dictionary of arrays that hold each subset array
    subsetArrays = {}

    # Groups by the subset classifier and then creates groups from that 
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
            subsetArrays[group] = multinet.array.loc[interStructs, :, :]

        else:
            print(f"No structures in MultiNetwork for group: {group}")

    return subsetArrays