# %%
from my_functions import *

# %%
file= uproot.open("/home/ruben/Escritorio/BachelorThesisRuben/Data/Muon_trees.root") #opening the Root file with Uproot 

MuonTree_ZeroBias=file["MuonTree_ZeroBias;1"]

#For now I'm restricting the dataset to 3000
eta=MuonTree_ZeroBias.arrays(["LVL1Muon_eta"], entry_start=0, entry_stop=3000)
phi=MuonTree_ZeroBias.arrays(["LVL1Muon_phi"], entry_start=0, entry_stop=3000)

# %%
a=all_events_jTower_dr(eta,phi,MuonTree_ZeroBias,MuonTree_ZeroBias,"jTower_eta","jTower_phi")

a


