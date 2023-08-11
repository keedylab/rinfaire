from IndividualNetwork import IndividualNetwork
from Structure import Structure

struct = Structure('/data/araju/ptpsinpdb/FinalRefine_qFitV4_V2/6B8Z/6B8Z_qFit.pdb')
net = IndividualNetwork(struct)
net.populateNetwork()

print(struct.model, struct.name)
print(net.structure, net.name, net.network)