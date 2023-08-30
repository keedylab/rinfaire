import gemmi
import os

class Structure:
    def __init__ (self, pathname, args):
        
        self.setModel(pathname)
        self.setName(pathname)
        self.setSequence()
        self.setChains()

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

    def setChains (self):
        
        # Creates list of chains in structure
        self.chains = []

        # Iterates over all the chains in the model, append to this list
        for n_ch, chain in enumerate(self.model[0]):
            self.chains.append(Chain(chain, self))

    def setSequence (self):
        self.sequence = {}

        # Iterates over all chains and residues, asserts that the residue is not a HETATM
        for n_ch, chain in enumerate(self.model[0]):
            for n_res, res in enumerate(chain):

                perResiAtomCouter = 0

                if res.het_flag == 'A':
                    for n_atom, atom in enumerate(res):

                        perResiAtomCouter += 1
                    
                    # Appends to the sequence dict the:
                    # Residue number as key and amino acid as value
                    self.sequence[res.seqid.num] = {}
                    self.sequence[res.seqid.num]['name'] = res.name
                    self.sequence[res.seqid.num]['atomcount'] = perResiAtomCouter

                    #print(res.seqid.num, res.name, perResiAtomCouter)

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
    
    def getChains (self):
        return self.chains
    
    # def getPeptideChains (self):
        
    #     # Interates
    #     for chain in self.chains:
    #         for n_res, res in enumerate(chain):
                
    def getSequence (self):
        return self.sequence
    
    def getFirstResi (self):

        # Gets keys of self.sequence and creates a list with them
        # Since dictionaries are stored by insertion order and the first item inserted is the first residue,
        # First element ([0]) is the first residue
        startingResi = list(self.sequence.keys())[0]
        return startingResi

class Chain:

    def __init__ (self, chain, structure):
        
        self.chain = chain
        self.structure = structure
