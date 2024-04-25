import xarray as xr
import numpy as np
import pandas as pd
from Bio import SeqIO
import pickle
import logging

class MultiNetwork:
    
    # Class constructor
    def __init__ (self, args=None, array=None, seqaln=None, metadata=None):

        """
        Function is the class constructor for MultiNetwork
        """
        
        # Sets args as class variable
        # Tests whether args are already provided or not
        if args is not None:
            self.args = args

        else:
            self.args = None

        # Sets the XArray multidimensional array associated with the MultiNetwork
        # Tests whether array is already provided or not
        if array is not None:
            self.array = array

        else:
            self.array = None

        # Sets the sequence alignment
        # Tests whether alignment is already provided or not
        if seqaln is not None:
            self.seqaln = seqaln

        # Condition if alignment file is provided in args
        elif args.alignmentFile is not None:
            self.setSeqAlignment(args.alignmentFile)
        
        else:
            self.seqaln = None

        # Sets the metadata class variable
        # Condition if metadata is provided as optional argument
        if metadata is not None:
            self.metadata = metadata
        
        # Condition if metadata file is provided in args
        elif args.metadata is not None:
            self.setMetaData(args.metadata)
        
        else:
            self.metadata = None

    # Set functions

    def setSeqAlignment (self, alignmentPath):
        # Uses the SeqIO tool from Biopython to parse the fasta file
        self.seqaln = {}
        self.size = 0
        
        for record in SeqIO.parse(alignmentPath, "fasta"):
            pdbid = record.id
            self.seqaln[pdbid] = record.seq
            
            # Makes sure that self.size is as long as the longest sequence (although they should be all the same) 
            if len(record.seq) > self.size:
                self.size = len(record.seq)
            
        # Increments by one due to how indexing starts at 0 but residue numbering starts at 1
        self.size += 1

    def setMetaData (self, csvFile):

        """
        Function that takes in a csv file of associated metadata and adds it as a pandas df to be a class variable of MultiNetwork
        """

        # Imports the csv and stores it as Pandas dataframe
        try:
            self.metadata = pd.read_csv(csvFile, header=0)
        except FileNotFoundError:
            print(f"File not found: {csvFile}")
            return
        except Exception as e:
            print(f"Error reading CSV file: {e}")
            return
           
        # Converts all entries in pandas table into a string (for now)
        self.metadata = self.metadata.astype(str)

        # Iterates over each column
        for metadataColumn in self.metadata.columns:

            # Splits values delimited as ; into a list
            kwargs = {metadataColumn : self.metadata[metadataColumn].str.split('; ')}
            self.metadata = self.metadata.assign(**kwargs)

    def oneToAll (self, seqID, sequenceList, seqResidue):
        
        mainCount = 0
        seqIndex = 0
     
        # Iterates over each residue in the sequence of the structure being queried
        for i in self.seqaln[seqID]:

            # Increments the main count that counts the total number of positions it has passed on the alignment
            mainCount = mainCount + 1

            # If there is a residue present (not just a - used as a placeholder)
            # Then increment the seq index that counts the total number of actual residue positions that have been passed
            if i != '-':
                seqIndex = seqIndex + 1

            # Condition when it reaches the residue in question 
            # (since seq index represents the position in the individual sequence and not the full alignment)
            if sequenceList[seqIndex] == int(seqResidue) and i != '-':
                
                logging.debug(f'Final Main count: {mainCount} Seq count: {sequenceList[seqIndex]} Resi: {i}')

                # Returns the main count which is the position on the full alignment
                return(mainCount)

            logging.debug(f'Main count: {mainCount} Seq count: {sequenceList[seqIndex]} Resi: {i}')

    def add (self, inputAdjacencyDict, inputStructName, inputSequenceList):

        # Loops through all pairings                
        for firstResi in inputAdjacencyDict:
            for secondResi in inputAdjacencyDict[firstResi]:
                
                #print(f"Before conversion: {firstResi} {secondResi}")

                # Uses the OneToAll conversion function to map individual resi # to sequence alignment #
                updatedFirstResi = self.oneToAll(inputStructName, inputSequenceList, firstResi)
                updatedSecondResi = self.oneToAll(inputStructName, inputSequenceList, secondResi)

                #print(f"After conversion: {updatedFirstResi} {updatedSecondResi}")

                logging.debug(f'Adding network pair: ({firstResi},{secondResi}) as ({updatedFirstResi},{updatedSecondResi})')
                # Adds pairing to the array along with the weight
                self.array.loc[inputStructName, updatedFirstResi, updatedSecondResi] = inputAdjacencyDict[firstResi][secondResi]['weight']

    # TODO: Update unit test to make sure this function works
    def normalizeStruct (self):

        # Creates a vector of maximum values across the first and second residue
        # Essentially a maximum value for each network
        maxValues = self.array.max(dim=['firstResi','secondResi'])

        # Then divides each network by the corresponding value in the vector of max values
        # Scales from 0 - 10
        self.array = (self.array / maxValues[:]) * 10

    def scaleMultiNet (self):

        # Gets maximum value across all dimensions
        maxValue = self.array.max(dim=['network','firstResi','secondResi']).item()

        # Then divides each network by max value
        # Scales to a set value (default 0 to 20)
        self.array = (self.array / maxValue) * self.args.multinet_scale
    
    # TODO: Update unit test to make sure this function works
    def addNetworks (self, networkList):

        # Creates list to represent all the structures
        structList = []
        for net in networkList:
            structList.append(net.struct.getName())

        print(structList)

        # Creates new blank 3D array that is of size (number of structures x length of seq x length of seq)
        self.array = xr.DataArray(
            0.0, 
            coords=dict(network=structList, firstResi=range(self.size), secondResi=range(self.size)), 
            dims=("network", "firstResi", "secondResi")
        )

        # Iterates over the network list and adds values from each adjacency matrix
        for net in networkList:
            logging.info(f'Adding network: {net.struct.getName()}')
            self.add(
                net.convertToAdjacency(), 
                net.struct.getName(), 
                net.struct.getSequenceList()
            )

        # Normalization of edge weights relative to the whole structure being added
        # All values are scaled from (0) to the max edge weight present (10)
        if self.args.no_norm_struct == False:
            self.normalizeStruct()
            logging.info(f'Normalized all individual structures')

        # Scales the entire MultiNetwork to be between values 0 and 10
        if self.args.scale_multinet == True:
            self.scaleMultiNet()
            logging.info(f'Scaled the MultiNetwork to values between 0 and 10')

        # Replaces NaN values with zeroes
        self.array = self.array.fillna(0)

        logging.info(f'Finished adding networks to MultiNetwork object')
        print(self.array)

    def getInfo (self):

        self.getInfo_Edges()
        self.getInfo_Structs()

    def getInfo_Edges (self):

        import matplotlib.pyplot as plt

        # Gets array into 1D list-like form, and then removes instances where's there's no value (0)
        stackedArray = self.array.stack(allDims=[...])
        stackedArray = stackedArray.where(stackedArray > 0, drop=True)

        minValue = stackedArray.min().values
        maxValue = stackedArray.max().values

        print(f'Minimum value in array is: {minValue}')
        print(f'Maximum value in array is: {maxValue}')

        # Plots histogram of edge distribution across the network
        outputInfoName = f'{self.args.output}MultiNetwork_InfoEdges'
        plt.figure(figsize=(10,10))
        xr.plot.hist(stackedArray, bins=range(0, int(maxValue) + 1, 1), range=(0, maxValue))
        plt.savefig(outputInfoName + '.png')
        plt.clf()

    def getInfo_Structs (self):

        import matplotlib.pyplot as plt

        # Sums across both residue axes to get the sum value for each network
        summedArray = self.array.sum(dim=['firstResi','secondResi'])

        minValue = summedArray.min().values
        maxValue = summedArray.max().values

        print(f'Lowest weight network is: {summedArray.idxmin().values} with value of: {minValue}')
        print(f'Highest weight network is: {summedArray.idxmax().values} with value of: {maxValue}')

        # Plots histogram of network weight distribution
        outputInfoName = f'{self.args.output}MultiNetwork_InfoStructs_Hist'
        plt.figure(figsize=(10,10))
        xr.plot.hist(summedArray, bins=range(0, int(maxValue) + 10, 10), range=(0, maxValue))
        plt.savefig(outputInfoName + '.png')
        plt.clf()

        # Plots network weights by descending order
        plt.figure(figsize=(50,50))
        networkSeries = summedArray.to_series().sort_values(ascending=False)
        networkSeries.plot.bar()

        outputInfoName = f'{self.args.output}MultiNetwork_InfoStructs_Bar'
        plt.savefig(outputInfoName + '.png')
        plt.clf()

    def exportPickle (self):
        
        # Creates new pickle (.pkl) file and then dumps the entire class object into the pickle file
        with open(f'{self.args.output}MultiNetwork.pkl', 'wb') as pickleFile:
            pickle.dump(self, pickleFile)