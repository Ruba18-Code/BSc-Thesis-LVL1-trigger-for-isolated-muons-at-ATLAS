# %%
from my_functions import *

# %%
file= uproot.open("/home/ruben/Escritorio/BachelorThesisRuben/Data/Muon_trees.root") #opening the Root file with Uproot 

MuonTree_ZeroBias=file["MuonTree_ZeroBias;1"]

# %%
eta=MuonTree_ZeroBias.arrays(["LVL1Muon_eta"], entry_start=2000, entry_stop=2100)
phi=MuonTree_ZeroBias.arrays(["LVL1Muon_phi"], entry_start=2000, entry_stop=2100)

# %%
a=all_events_jTower_dr(eta,phi,MuonTree_ZeroBias,MuonTree_ZeroBias,"jTower_eta","jTower_phi")

a

# %%
print(isAll_below_threshold(a,0.4))


