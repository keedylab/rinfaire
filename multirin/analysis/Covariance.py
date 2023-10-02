import xarray as xr
import networkx as nx
from pyvis.network import Network
import logging
import pickle

class Covariance:

    def __init__ (self, args):
        self.args = args

    def readPickle (self):
        
        # Opens pickle file
        with open(self.args.filename, 'rb') as pickleFile:
            self.multinet = pickle.load(pickleFile)

    def flatten (self):
        self.flattenMulti = self.multinet.array.stack(resiPair=("firstResi", "secondResi"))

    def exportPickle (self):

        # Creates new pickle (.pkl) file and then dumps the entire class object into the pickle file
        with open(f'{self.args.outputname}.pkl', 'wb') as pickleFile:
            pickle.dump(self, pickleFile)
