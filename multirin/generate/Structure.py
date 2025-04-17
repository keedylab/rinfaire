import gemmi
import os
import logging

class Structure:
    def __init__ (self, pathname, args):
        
        self.setModel(pathname)
        self.setName(pathname)
        self.setSequence()
        self.setSequenceList()

        self.args = args
        
        self.spacegroup = None

    # Set functions

    def setModel (self, pathname):
        # Reads structure through Gemmi, sets up entities
        # Then sets this to the model variable in the class
        st = gemmi.read_structure(pathname)
        st.setup_entities()
        self.model = st

    def setName (self, pathname):
        # Sets the name of the structure to be the filename
        filename = os.path.basename(pathname)
        filename = filename[:filename.rindex('.')]
        self.name = filename

    def setSequence (self):
        self.sequence = {}

        # Iterates over all chains and residues, asserts that the residue is not a HETATM
        for n_ch, chain in enumerate(self.model[0]):
            for n_res, res in enumerate(chain):

                perResiAtomCounter = 0

                if res.het_flag == 'A':
                    for n_atom, atom in enumerate(res):

                        perResiAtomCounter += 1
                    
                    # Appends to the sequence dict the:
                    # Residue number as key and amino acid as value
                    self.sequence[res.seqid.num] = {}
                    self.sequence[res.seqid.num]['name'] = res.name
                    self.sequence[res.seqid.num]['atomcount'] = perResiAtomCounter

                    logging.debug(f'Number:{res.seqid.num} ResName:{res.name} AtomCount:{perResiAtomCounter}')

    def setSequenceList (self):
        
        # Creates list with one element (so that index 0 will have a value)
        self.sequenceList = [0]

        # Then adds list of sequence keys to this list
        self.sequenceList = self.sequenceList + sorted(list(self.sequence.keys()))

    # def setAttributes (self):
    #     cifBlock = self.model.make_mmcif_headers()

    #     print(cifBlock.get_mmcif_category('_diffrn'))
    #     print(cifBlock.get_mmcif_category('_symmetry')['space_group_name_H-M'][0])
    #     print(cifBlock.get_mmcif_category('_entity'))

    # Get functions

    def getModel (self):
        return self.model
    
    def getName (self):
        return self.name

    def getSequence (self):
        return self.sequence

    def getSequenceList (self):
        return self.sequenceList
    
    def getFirstResi (self):

        # Gets keys of self.sequence and creates a list with them
        # Since dictionaries are stored by insertion order and the first item inserted is the first residue,
        # First element ([0]) is the first residue
        startingResi = self.sequenceList[1]
        return startingResi

        # startingResi = list(self.sequence.keys())[0]
        # return startingResi