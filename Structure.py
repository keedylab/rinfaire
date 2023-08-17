import gemmi
import os

class Structure:
    def __init__ (self, pathname):
        
        self.setModel(pathname)
        self.setName(pathname)
        self.setSequence()
        
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
        self.sequence = []

        # Iterates over all chains and residues, asserts that the residue is not a HETATM
        for n_ch, chain in enumerate(self.model[0]):
            for n_res, res in enumerate(chain):
                if res.het_flag == 'A':
                    
                    # Appends to the sequence list the:
                    # Residue number and amino acid as a tuple
                    self.sequence.append((res.seqid.num, res.name))

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
    
    def getFirstResi (self):
        return self.sequence[0][0]