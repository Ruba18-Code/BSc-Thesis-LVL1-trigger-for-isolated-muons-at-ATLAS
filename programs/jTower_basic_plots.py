# %%
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from my_functions import *

# %%
file= uproot.open("/home/ruben/Escritorio/BachelorThesisRuben/Data/Muon_trees.root") #opening the Root file with Uproot 

MuonTree_ZeroBias=file["MuonTree_ZeroBias;1"]

# %%
xmin=10**4
xmax=8*10**4
binsize=abs(xmin-xmax)/100
chunk_size=1000

xlabel=r"$e_T$ (MeV)"
ylabel=r"Counts"
title=r"$e_T$ vs counts - jTower data"

tree=MuonTree_ZeroBias
name="jTower_et_MeV"

_,_, = jTower_handler(tree,name,xmin,xmax,binsize,chunk_size,xlabel,ylabel,title,True)

# %%
xmin=-5
xmax=5
binsize=abs(xmin-xmax)/100
chunk_size=1000

xlabel=r"$\eta$"
ylabel=r"Counts"
title=r"$\eta$ vs counts - jTower data"

tree=MuonTree_ZeroBias
name="jTower_eta"

_,_, = jTower_handler(tree,name,xmin,xmax,binsize,chunk_size,xlabel,ylabel,title,True)

# %%
xmin=-4
xmax=4
binsize=abs(xmin-xmax)/100
chunk_size=1000

xlabel=r"$\phi$"
ylabel=r"Counts"
title=r"$\phi$ vs counts - jTower data"

tree=MuonTree_ZeroBias
name="jTower_phi"

_,_, = jTower_handler(tree,name,xmin,xmax,binsize,chunk_size,xlabel,ylabel,title,True)


