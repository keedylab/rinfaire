import gemmi
import os

class Structure:
    def __init__ (self, pathname):
        st = gemmi.read_structure(pathname)
        st.setup_entities()
        self.model = st
        
        filename = os.path.basename(pathname)
        filename = filename[:filename.rindex('.')]
        self.name = filename

        self.spacegroup = None

    

