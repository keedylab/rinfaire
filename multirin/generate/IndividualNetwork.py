#!/usr/bin/env python
# coding: utf-8

import gemmi
import networkx as nx
from pyvis.network import Network
import copy
import logging

class IndividualNetwork:
    def __init__ (self, Structure, args):
        
        self.struct = Structure
        self.network = nx.Graph() # Creates networkX graph object
        self.args = args
    
    def findAltConfAtoms (self):

        # Creates empty dictionary of atoms with alt-confs
        atomsWithAltConfsDict = {}

        # Iterates over all chains, residues, and atoms to only add atoms that have an altloc
        for n_ch, chain in enumerate(self.struct.model[0]):
            for n_res, res in enumerate(chain):
                
                # Asserts that the residue is part of the polymer and not a ligand or water (by ensuring that it's not a HETATM)
                if res.het_flag == 'A':
                    #print("Polymer: ", res.name,res,chain)
                    
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

        return atomsWithAltConfsDict

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
    
    def updateEdge (self, firstResi, secondResi, counterResiResi, backboneAtoms=None):

        # Conditional based on if normalizing the atom-atom pair numbers based on the residue type is toggled on/off
        # Normally normalization is: ON 
        if self.args.no_norm_resi == False:

            # Normalizes by taking each atom-atom pair added and dividing it by the total # of atoms
            # Equivalent to taking the full count of atom-atom pairs for this residue pair and dividing by total # of atoms

            # Condition if this is called from backbone
            if backboneAtoms is not None:
                # avgNumberAtomsInResidue = 15
                addValue = 10 / (len(backboneAtoms) * 4)
            else:
                addValue = 10 / (self.struct.sequence[firstResi]['atomcount'] + self.struct.sequence[secondResi]['atomcount'])
        else:
            addValue = 0.1

        # Tests if the connection it found is already present in the edge list
        # If so, it increments the edge weight by 0.1. If not, it creates a new edge
        if (firstResi,secondResi) in self.network.edges:
            #print('Adding to existing connection between: ', firstResi, secondResi)
            self.network[firstResi][secondResi]['weight'] = self.network[firstResi][secondResi]['weight'] + addValue
        else:
            #print('Creating new connection between: ', firstResi, secondResi)
            self.network.add_edge(firstResi, secondResi, weight=addValue) # Adds a new edge with a weight of 0.1
            counterResiResi += 1 #Increments the count of residue-residue connections by 1
        
        return counterResiResi
    
    def populateNetwork (self):
        atomsWithAltConfsDict = self.findAltConfAtoms()
        amideHOnlyList = self.flagAmideHydrogenOnlyResidues(atomsWithAltConfsDict)
        
        contactCutoffValue = 4 # Distance cutoff value
        tooFarCutoffValue = 25 # Minimum distance for two atoms and therefore residues to be considered as too far from each other

        counterResiResi = 0

        #List of all backbone atoms in the protein structure
        backboneAtoms = ['N','CA','C','O'] # 'H','HA','HA2','HA3'

        for firstResi in atomsWithAltConfsDict:
            for secondResi in atomsWithAltConfsDict:

                # Condition to remove any cases of residues with amide H alt confs having connections with adjacent residues on the backbone
                if ((firstResi in amideHOnlyList) or (secondResi in amideHOnlyList)) and (firstResi + 1 == secondResi):
                    #print("These residues' connections:", firstResi, secondResi, "are not being searched because they are adjacent and one has amide H's")
                    continue

                # Condition that satisfies both the fact that the first and second residues cannot be equal to each other
                # And that we can prune duplicate connections by only looking at connection i,j and not j,i
                if firstResi < secondResi:

                    tooFarFlag = False

                    # Finds atom names for all atoms in first and second residues
                    atomsWithAltConfsNamesFirstResi = []
                    for atom in atomsWithAltConfsDict[firstResi]:
                        atomsWithAltConfsNamesFirstResi.append(atom.name)

                    atomsWithAltConfsNamesSecondResi = []
                    for atom in atomsWithAltConfsDict[secondResi]:
                        atomsWithAltConfsNamesSecondResi.append(atom.name)

                    # First condition to test if the two residues are adjacent residues that have coupled backbone alt-confs
                    if (firstResi + 1 == secondResi) and (set(backboneAtoms).issubset(set(atomsWithAltConfsNamesSecondResi)) and set(backboneAtoms).issubset(set(atomsWithAltConfsNamesFirstResi))):

                        print(atomsWithAltConfsNamesSecondResi)
                        backboneAtomsFirstResi = [i in backboneAtoms for i in atomsWithAltConfsNamesFirstResi]
                        filtered = list(filter(lambda x: backboneAtomsFirstResi, atomsWithAltConfsNamesFirstResi))
                        print(filtered)

                        print(firstResi, secondResi, self.struct.name)

                        for firstAtom in backboneAtoms:
                            for secondAtom in backboneAtoms:
                                #print("Found backbone connection between: ", firstResi, secondResi, firstAtom, secondAtom)
                                counterResiResi = self.updateEdge(firstResi, secondResi, counterResiResi,backboneAtoms=backboneAtoms)
                    
                    # All other cases when not backbone alt confs
                    else:
                        for firstAtom in atomsWithAltConfsDict[firstResi]:
                            for secondAtom in atomsWithAltConfsDict[secondResi]:

                                # If not, then calculate the distance between the two atoms
                                # Calculates distance between the two atoms using Gemmi dist() function
                                distance = firstAtom.pos.dist(secondAtom.pos)  

                                # Condition that checks if the distance between atoms is greater than the threshold
                                # If it is, then they are considered to be too far to even both checking the rest of the residue
                                if distance > tooFarCutoffValue:
                                    tooFarFlag = True
                                    #print("These residues are too far to warrant a full atom-atom search: ", firstResi, secondResi, "checking next connection")
                                    break
                                
                                # Asks if the distance calculated between each atom pair is less than the maximum atomic distance the user specifies
                                elif distance < contactCutoffValue:
                                    #print("Found distance connection between: ", firstResi, secondResi, firstAtom, secondAtom)
                                    counterResiResi = self.updateEdge(firstResi, secondResi, counterResiResi)

                            # If the tooFarFlag is triggered then it continues to break this loop to prevent it from searching any atom-atom contacts...
                            # ...between this pair and move on to the next pair of residues
                            if tooFarFlag == True:
                                break
    
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



