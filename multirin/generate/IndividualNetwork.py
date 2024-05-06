#!/usr/bin/env python
# coding: utf-8

import gemmi
import networkx as nx
from pyvis.network import Network
import copy
import logging
import numpy as np

class IndividualNetwork:
    def __init__ (self, Structure, args):
        
        self.struct = Structure
        self.network = nx.Graph() # Creates networkX graph object
        self.args = args
    
    def findAltConfAtoms (self):

        print(f'Starting structure: {self.struct.name}')

        # Creates empty dictionary of atoms with alt-confs
        atomsWithAltConfsDict = {}

        # Creates dictionary to map residue object in Gemmi with the residue number
        resToSeqPositionMap = {}

        # Iterates over all chains, residues, and atoms to only add atoms that have an altloc
        for n_ch, chain in enumerate(self.struct.model[0]):
            for n_res, res in enumerate(chain):
                
                # Asserts that the residue is part of the polymer and not a ligand or water (by ensuring that it's not a HETATM)
                if res.het_flag == 'A':
                    #print("Polymer: ", res.name,res,chain)

                    resToSeqPositionMap[res.seqid.num] = res
                    
                    # Iterates over each atom in the residue
                    for n_atom, atom in enumerate(res):
                
                        # Checks if the atom has an altloc label
                        if atom.has_altloc():

                            #If it does then append a tuple of (atom,residue) to the list of atoms with alt confs
                            if res.seqid.num in atomsWithAltConfsDict.keys():
                                #print('Resi present: ', res)
                                atomsWithAltConfsDict[res.seqid.num].append(atom)
                                
                            else:
                                #print('New resi: ', res)
                                atomsWithAltConfsDict[res.seqid.num] = []
                                atomsWithAltConfsDict[res.seqid.num].append(atom)

        return atomsWithAltConfsDict, resToSeqPositionMap

    def flagAmideHydrogenOnlyResidues (self, atomsWithAltConfsDict):
        amideHOnlyList = []

        for altResi in atomsWithAltConfsDict.copy():

            listAtoms = []
            for altAtoms in atomsWithAltConfsDict[altResi]:
                listAtoms.append(altAtoms.name)

            # Uses a set to ensure that all of the elements present in the list are all Hydrogens
            # Since the number of Hydrogens could vary based on how many altlocs there are
            if set(listAtoms) == set('H'):
                amideHOnlyList.append(altResi)
                #print("Adding to amide Hydrogen only list: ", altResi, " because it only has amide hydrogen alt-confs: ", listAtoms)

        return amideHOnlyList
    
    def populateNetwork (self):
        atomsWithAltConfsDict, resToSeqPositionMap = self.findAltConfAtoms()
        amideHOnlyList = self.flagAmideHydrogenOnlyResidues(atomsWithAltConfsDict)
        
        # Creates a directed graph representing the path of search for backbone atoms in coupled backbone instances
        backboneGraph = nx.DiGraph({'N_2': ['C_1','H_2','CA_2'], 
                                    'CA_2': ['C_2','HA_2'],
                                    'C_2': ['O_2'],
                                    'C_1': ['O_1','CA_1'],
                                    'CA_1': ['HA_1','N_1'],
                                    'N_1': ['H_1']
                                    })

        # Creates lists for tracking the weights of all the following types of connections
        self.weightsRecord = {'adjResi': {'total': [], 'BB_BB': [], 'SC_BB': [], 'SC_SC': []},
                              'nonAdjResi': {'total': []}}
        
        self.distancesRecord = {'adjResi': {'total': [], 'SC_BB': [], 'SC_SC': []},
                              'nonAdjResi': {'total': []}}

        # Iterates over all pairs of residues that have alt-confs
        for firstResi in atomsWithAltConfsDict:
            for secondResi in atomsWithAltConfsDict:

                # Condition that satisfies both the fact that the first and second residues cannot be equal to each other
                # And that we can prune duplicate connections by only looking at connection i,j and not j,i
                if firstResi < secondResi:

                    # Sets total connections
                    totalConnections = 0

                    # Sets the normalization factor
                    normalizationFactor = (self.struct.sequence[firstResi]['atomcount'] + self.struct.sequence[secondResi]['atomcount']) / 10

                    # Condition to remove any cases of residues with amide H alt confs having connections with adjacent residues on the backbone
                    if ((firstResi in amideHOnlyList) or (secondResi in amideHOnlyList)) and (firstResi + 1 == secondResi):
                        #print("These residues' connections:", firstResi, secondResi, "are not being searched because they are adjacent and one has amide H's")
                        continue

                    # Find backbone vs sidechain atoms
                    backboneGraphCopy = copy.deepcopy(backboneGraph)

                    # Conditions to either add or remove edges if either first or second residue is a Proline or Glycine
                    if self.struct.sequence[firstResi]['name'] == 'PRO':
                        backboneGraphCopy.remove_edge('N_1','H_1')

                    if self.struct.sequence[secondResi]['name'] == 'PRO':
                        backboneGraphCopy.remove_edge('N_2','H_2')

                    if self.struct.sequence[firstResi]['name'] == 'GLY':
                        backboneGraphCopy.remove_edge('CA_1','HA_1')
                        backboneGraphCopy.add_edge('CA_1','HA2_1')
                        backboneGraphCopy.add_edge('CA_1','HA3_1')

                    if self.struct.sequence[secondResi]['name'] == 'GLY':
                        backboneGraphCopy.remove_edge('CA_2','HA_2')
                        backboneGraphCopy.add_edge('CA_2','HA2_2')
                        backboneGraphCopy.add_edge('CA_2','HA3_2')

                    # Creates lists of backbone vs sidechain atoms for both residues, also creates the same list but just for the atom names
                        
                    backboneAtomsFirstResi, sidechainAtomsFirstResi, backboneAtomsSecondResi, sidechainAtomsSecondResi = [], [], [], []
                    backboneAtomsFirstResiNames, sidechainAtomsFirstResiNames, backboneAtomsSecondResiNames, sidechainAtomsSecondResiNames = [], [], [], []
                    backboneAtomsFirstResiGroupedByAltLoc, sidechainAtomsFirstResiGroupedByAltLoc, backboneAtomsSecondResiGroupedByAltLoc, sidechainAtomsSecondResiGroupedByAltLoc = {}, {}, {}, {}
                    backboneAtomsFirstResiNamesGroupedByAltLoc, sidechainAtomsFirstResiNamesGroupedByAltLoc, backboneAtomsSecondResiNamesGroupedByAltLoc, sidechainAtomsSecondResiNamesGroupedByAltLoc = {}, {}, {}, {}

                    # Iterates over each atom in the first residue, checks if atom name is in backbone graph, if so then adds to backbone atoms and if not adds to sidechain atoms
                    for atom in atomsWithAltConfsDict[firstResi]:
                        if f'{atom.name}_1' in backboneGraphCopy.nodes:
                            backboneAtomsFirstResi.append(atom)
                            backboneAtomsFirstResiNames.append(f'{atom.name}_1')

                            if atom.altloc not in backboneAtomsFirstResiGroupedByAltLoc:
                                backboneAtomsFirstResiGroupedByAltLoc[atom.altloc] = []
                                backboneAtomsFirstResiNamesGroupedByAltLoc[atom.altloc] = []
                            backboneAtomsFirstResiGroupedByAltLoc[atom.altloc].append(atom)
                            backboneAtomsFirstResiNamesGroupedByAltLoc[atom.altloc].append(f'{atom.name}_1')
                            
                        else:
                            sidechainAtomsFirstResi.append(atom)
                            sidechainAtomsFirstResiNames.append(f'{atom.name}_1')

                            if atom.altloc not in sidechainAtomsFirstResiGroupedByAltLoc:
                                sidechainAtomsFirstResiGroupedByAltLoc[atom.altloc] = []
                                sidechainAtomsFirstResiNamesGroupedByAltLoc[atom.altloc] = []
                            sidechainAtomsFirstResiGroupedByAltLoc[atom.altloc].append(atom)
                            sidechainAtomsFirstResiNamesGroupedByAltLoc[atom.altloc].append(f'{atom.name}_1')

                    # Does the same but for the second residue
                    for atom in atomsWithAltConfsDict[secondResi]:
                        if f'{atom.name}_2' in backboneGraphCopy.nodes:
                            backboneAtomsSecondResi.append(atom)
                            backboneAtomsSecondResiNames.append(f'{atom.name}_2')

                            if atom.altloc not in backboneAtomsSecondResiGroupedByAltLoc:
                                backboneAtomsSecondResiGroupedByAltLoc[atom.altloc] = []
                                backboneAtomsSecondResiNamesGroupedByAltLoc[atom.altloc] = []
                            backboneAtomsSecondResiGroupedByAltLoc[atom.altloc].append(atom)
                            backboneAtomsSecondResiNamesGroupedByAltLoc[atom.altloc].append(f'{atom.name}_2')

                        else:
                            sidechainAtomsSecondResi.append(atom)
                            sidechainAtomsSecondResiNames.append(f'{atom.name}_2')

                            if atom.altloc not in sidechainAtomsSecondResiGroupedByAltLoc:
                                sidechainAtomsSecondResiGroupedByAltLoc[atom.altloc] = []
                                sidechainAtomsSecondResiNamesGroupedByAltLoc[atom.altloc] = []
                            sidechainAtomsSecondResiGroupedByAltLoc[atom.altloc].append(atom)
                            sidechainAtomsSecondResiNamesGroupedByAltLoc[atom.altloc].append(f'{atom.name}_2')

                    # Remove backbone alt confs from calculation if specified (by default we do look at backbones and this is FALSE)
                    if self.args.only_sidechain == True:

                        # Finds connections only between sidechain atoms in first residue and second residue
                        totalConnections = self.findConnections(sidechainAtomsFirstResi, sidechainAtomsSecondResi)

                    # Includes backbone alt confs (normal running scenario)
                    else:

                        # Case when these are adjacent residues
                        if (firstResi + 1 == secondResi):

                            # STEP 1: Check backbone atoms for connections with other backbone atoms in continuous stretch of alt conf atoms
                            
                            # Checks if there are alt confs across the peptide bond (first residue has a carbonyl C alt-conf, second residue as a backbone amide N alt-conf)
                            # If they don't then there are no backbone connections
                            if (('C_1' in backboneAtomsFirstResiNames) and ('N_2' in backboneAtomsSecondResiNames)) == False:
                                BB_BB_Connections = 0

                            # If they do both have alt confs, use recursive search to find number of connections
                            else:
                                BB_BB_Connections = 0

                                # Finds common alt-loc labels between the two residues
                                altLocsInFirstResi = set(backboneAtomsFirstResiGroupedByAltLoc.keys())
                                altLocsInSecondResi = set(backboneAtomsSecondResiGroupedByAltLoc.keys())
                                altLocsInBoth = altLocsInFirstResi.intersection(altLocsInSecondResi)

                                # For each common alt-loc, run the recursive backbone search to find backbone connections across the two residues
                                for altloc in altLocsInBoth:
                                    BB_BB_Connections += self.findBackboneConnections(backboneAtomsFirstResiNamesGroupedByAltLoc[altloc], backboneAtomsSecondResiNamesGroupedByAltLoc[altloc], backboneGraphCopy)
                                
                                # Normalizes the number of connections by taking the connections and dividing by the number of alt conf atoms in both first and second residue (also scales by factor of 10 for visualization)

                                # This allows us to normalize for both:
                                # - Larger residues (eg. Trp) having more atoms that could bias the calculation to highlight larger residues
                                # - Cases where there are a lot of alt confs for a given residue (eg. alt A,B,C,D) which would inflate the total number of atom-atom connections
                                
                                # Normally normalization is: ON, therefore this is FALSE
                                if self.args.no_norm_resi == False:
                                    # backboneConnections = (backboneConnections / (len(set(backboneAtomsFirstResiNames)) + len(set(backboneAtomsSecondResiNames)))) * 10
                                    BB_BB_Connections = BB_BB_Connections / normalizationFactor

                            # STEP 2: Check sidechain atoms for connections with other sidechain atoms in adjacent residue

                            # If there are sidechain alt confs present
                            if (len(sidechainAtomsFirstResi) + len(sidechainAtomsSecondResi)) > 0:
                                SC_SC_Connections, SC_SC_distancesRecord = self.findConnections(sidechainAtomsFirstResi,sidechainAtomsSecondResi)
                                SC_BB_Connections, SC_BB_distancesRecord = self.findConnections(sidechainAtomsFirstResi,backboneAtomsSecondResi) #excludeAtoms=['CB']
                                BB_SC_Connections, BB_SC_distancesRecord = self.findConnections(backboneAtomsFirstResi,sidechainAtomsSecondResi) #excludeAtoms=['CB']
                                
                                # Updates distance records
                                self.distancesRecord['adjResi']['SC_SC'] += SC_SC_distancesRecord
                                self.distancesRecord['adjResi']['SC_BB'] += SC_BB_distancesRecord + BB_SC_distancesRecord
                                self.distancesRecord['adjResi']['total'] += SC_SC_distancesRecord + SC_BB_distancesRecord + BB_SC_distancesRecord

                                # Same normalization as above
                                if self.args.no_norm_resi == False:
                                    SC_SC_Connections = SC_SC_Connections / normalizationFactor
                                    SC_BB_Connections = SC_BB_Connections / normalizationFactor
                                    BB_SC_Connections = BB_SC_Connections / normalizationFactor

                            if BB_BB_Connections != 0:
                                self.weightsRecord['adjResi']['BB_BB'].append(BB_BB_Connections)
                            if (SC_BB_Connections + BB_SC_Connections) != 0:    
                                self.weightsRecord['adjResi']['SC_BB'].append(SC_BB_Connections + BB_SC_Connections) # Since BB-SC and SC-BB are the same case just swapped in residue order
                            if SC_SC_Connections != 0:
                                self.weightsRecord['adjResi']['SC_SC'].append(SC_SC_Connections)

                            # Adds together the backbone and sidechain connections
                            totalConnections = BB_BB_Connections + BB_SC_Connections + SC_BB_Connections + SC_SC_Connections

                            # # Same normalization as above
                            # if self.args.no_norm_resi == False:
                            #     totalConnections = totalConnections / normalizationFactor

                            # Adds connection in the network between the two residues, with the weight being the total atom-atom connections
                            if totalConnections != 0:
                                self.network.add_edge(firstResi, secondResi, weight=totalConnections)

                                # Appends these weights to overall list of weights for tracking and subsequent plotting
                                # self.weightsRecord['adjResi']['BB_BB'].append(BB_BB_Connections)
                                # self.weightsRecord['adjResi']['SC_BB'].append(SC_BB_Connections + BB_SC_Connections) # Since BB-SC and SC-BB are the same case just swapped in residue order
                                # self.weightsRecord['adjResi']['SC_SC'].append(SC_SC_Connections)
                                self.weightsRecord['adjResi']['total'].append(totalConnections)

                        # Case when they are not adjacent residues
                        else:

                            # Finds connections between all atoms regardless of whether it's a backbone or sidechain atom 
                            totalConnections, nonAdj_distanceRecord = self.findConnections(atomsWithAltConfsDict[firstResi],atomsWithAltConfsDict[secondResi])

                            self.distancesRecord['nonAdjResi']['total'] += nonAdj_distanceRecord

                            # Same normalization as above
                            if self.args.no_norm_resi == False:
                                totalConnections = totalConnections / normalizationFactor

                            # Adds connection in the network between the two residues, with the weight being the total atom-atom connections
                            if totalConnections != 0:
                                self.network.add_edge(firstResi, secondResi, weight=totalConnections)
                            
                                # Appends this weight to tracking list
                                self.weightsRecord['nonAdjResi']['total'].append(totalConnections)

        # # Uses list comprehension to remove null values in these tracking/bookkeeping lists
        # self.adjResiBackboneWeights = [i for i in adjResiBackboneWeights if i != 0]
        # self.adjResiSidechainWeights = [i for i in adjResiSidechainWeights if i != 0]
        # self.adjResiTotalWeights = [i for i in adjResiTotalWeights if i != 0]
        # self.allAtomWeights = [i for i in allAtomWeights if i != 0]

        print(f"Average Adjacent Residue Total Weight for {self.struct.name}: {np.average(self.weightsRecord['adjResi']['total'])}")
        print(f"Average Non-Adjacent Residue Total Weight for {self.struct.name}: {np.average(self.weightsRecord['nonAdjResi']['total'])}")
    
    def findConnections (self, firstResiAltConfAtoms, secondResiAltConfAtoms, minDist=0, maxDist=4, tooFarDist=25, excludeAtoms=[]):

        """
        Function that finds distance connections between two residue's alt conf atoms.

        Inputs:
        - firstResiAltConfAtoms: List of alt conf atoms in first residue
        - secondResiAltConfAtoms: Same but for second residue
        - minDist: Minimum distance cutoff value (A)
        - maxDist: Maximum distance cutoff value (A)
        - tooFarDist: Minimum distance (A) for two atoms and therefore residues to be considered as too far from each other
        
        Outputs:
        - connections: Count of total number of connections
        - distancesRecord: List of the distances for each connection
        """

        connections = 0

        # List to keep track of distances
        distancesRecord = []

        for firstAtom in firstResiAltConfAtoms:
            for secondAtom in secondResiAltConfAtoms:

                # If not, then calculate the distance between the two atoms
                # Calculates distance between the two atoms using Gemmi dist() function
                distance = firstAtom.pos.dist(secondAtom.pos)  

                # Condition that checks if the distance between atoms is greater than the threshold
                # If it is, then they are considered to be too far to even both checking the rest of the residue
                if distance > tooFarDist:
                    # print("These residues are too far to warrant a full atom-atom search: ", firstResi, secondResi, "checking next connection")
                    return 0, []
                
                # Asks if the distance calculated between each atom pair is less than the maximum atomic distance the user specifies
                elif (distance < maxDist) and (distance > minDist):
                    if (firstAtom.name not in excludeAtoms) and (secondAtom.name not in excludeAtoms):
                        
                        connections += 1
                        distancesRecord.append(distance)

                    else:
                        print(firstAtom.name,secondAtom.name)
        
        return connections, distancesRecord

    def findBackboneConnections (self, firstResiAltConfAtomNames, secondResiAltConfAtomNames, backboneGraph, startingResidue='N_2', count=1):

        """
        """

        # Iterates over each successor node in the backbone graph from the current residue
        for successor in backboneGraph.successors(startingResidue):

            # If the successor has an alt-conf (present in either first residue or second residue alt conf names list)
            if (successor in firstResiAltConfAtomNames) or (successor in secondResiAltConfAtomNames):

                # Updates the count by 1
                newCount = count + 1

                # Recursively searches from the successor as the starting residue
                count = self.findBackboneConnections(firstResiAltConfAtomNames, secondResiAltConfAtomNames, backboneGraph, startingResidue=successor, count=newCount)

        return count
        
    
    def visualize (self):

        # Creates a deep copy of self.network
        # So any operations performed here do not affect the class object
        visNetwork = copy.deepcopy(self.network)

        # From VisArray function in postCONTACT script
        visNetwork.remove_nodes_from(list(nx.isolates(visNetwork)))
        # nx.convert_node_labels_to_integers(visNetwork)

        #Converts node labels into strings
        for i in visNetwork.nodes():
            visNetwork.nodes[i]['label'] = str(i)

        # Gets weights of network
        #widths = nx.get_edge_attributes(self.network, 'weight')

        # Creates network object
        # Notebook = true is for Jupyter Notebook (might need to remove later)
        nts = Network(notebook=True)

        # populates the nodes and edges data structures
        nts.from_nx(visNetwork)
        outputpath = f'{self.args.output}{self.struct.name}.html'
        nts.show(outputpath)

    def convertToAdjacency (self):
        return nx.to_dict_of_dicts(self.network)


# print("Total number of atom-atom connections is: ", counterAtomAtom)
# print("Total number of residue-residue connections is: ", counterResiResi)
# print(G.edges)
# print(G.edges.data('weight'))

# # Defining a separate function that looks for and labels residues that are adjacent to ligands

# # Creates empty dictionary of ligand atoms
# ligandAtomsDict = {}

# # Iterates over all chains, residues, and atoms to only add ligand atoms
# for n_ch, chain in enumerate(st[0]):
    
#     for n_res, res in enumerate(chain):
        
#         # Asserts that the atom is a HETATM (H) and that it is not a water -- should give us all ligand atoms
#         if res.het_flag == 'H' and (res.is_water() == False):
#             print("Ligand: ", res.name,res,chain)
            
#             # Iterates over each atom in the residue
#             for n_atom, atom in enumerate(res):

#                 #If it does then append a tuple of (atom,residue) to the list of atoms with alt confs
#                 if res.seqid.num in atomsWithAltConfsDict.keys():
#                     print('Ligand present: ', res)
#                     atomsWithAltConfsDict[res.seqid.num].append(atom)
                    
#                 else:
#                     print('New ligand: ', res)
#                     atomsWithAltConfsDict[res.seqid.num] = []
#                     atomsWithAltConfsDict[res.seqid.num].append(atom)

#                     # Checks if the atom has an altloc label
#                     # if atom.has_altloc():



