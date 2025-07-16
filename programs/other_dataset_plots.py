# %%
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from my_functions import*

# %%
file1= uproot.open("/home/ruben/Escritorio/BachelorThesisRuben/Data/MuonTree_ZB_0.root") #opening the Root file with Uproot 
file2=uproot.open("/home/ruben/Escritorio/BachelorThesisRuben/Data/Muon_trees.root")

# %%
MuonTree_ZeroBias=file1["MuonTree;1"]
MuonTree_Zmumu=file2["MuonTree_Zmumu"]
ZeroBias_pt=MuonTree_ZeroBias["muon_pt"].array()


dr_min=0.1
dr_max=0.3

nmin1=0
nmax1=6000
#Select quality 0 Z->mumu
Zmumu_pt=quality_selector_with_empty(MuonTree_Zmumu["muon_quality"].array(),MuonTree_Zmumu["muon_pt"].array(),0)[nmin1:nmax1]
Zmumu_eta=quality_selector_with_empty(MuonTree_Zmumu["muon_quality"].array(),MuonTree_Zmumu["muon_eta"].array(),0)[nmin1:nmax1]
Zmumu_phi=quality_selector_with_empty(MuonTree_Zmumu["muon_quality"].array(),MuonTree_Zmumu["muon_phi"].array(),0)[nmin1:nmax1]
#And select the Z peak pairs
Zmumu_pt, Zmumu_eta, Zmumu_phi = get_all_Z_peak_pairs(Zmumu_pt,Zmumu_eta,Zmumu_phi)
#Select the ZeroBias data with energy cut

nmin2=0
nmax2=100000
ZeroBias_eta=energy_cut_with_empty(MuonTree_ZeroBias["muon_pt"].array(), MuonTree_ZeroBias["muon_eta"].array())[nmin2:nmax2]
ZeroBias_phi=energy_cut_with_empty(MuonTree_ZeroBias["muon_pt"].array(), MuonTree_ZeroBias["muon_phi"].array())[nmin2:nmax2]
ZeroBias_pt=energy_cut_with_empty(MuonTree_ZeroBias["muon_pt"].array(), MuonTree_ZeroBias["muon_pt"].array())
print(len(ZeroBias_pt[ak.num(ZeroBias_pt) > 0]))

# %%
coolplot([Zmumu_pt, ZeroBias_pt], np.linspace(0,2*10**5, 50))

# %%
res1=muon_isolation_all_events(MuonTree_ZeroBias, ZeroBias_eta, ZeroBias_phi, dr_min, dr_max, [nmin2,nmax2])
res2=muon_isolation_all_events(MuonTree_Zmumu, Zmumu_eta, Zmumu_phi, dr_min, dr_max, [nmin1,nmax1])

# %%
coolplot([res1,res2], np.linspace(0,10000,50))

# %%
ratio1=ak.flatten(res1)/ak.flatten(ZeroBias_pt)
ratio2=ak.flatten(res2)/ak.flatten(Zmumu_pt)

coolplot([ratio1,ratio2], np.linspace(0,0.4,50))

# %%
dr_min=[0.1]
dr_max=[0.3]

plot_ROC_curve(MuonTree_Zmumu,MuonTree_ZeroBias,Zmumu_pt,Zmumu_eta,Zmumu_phi,ZeroBias_pt,ZeroBias_eta,ZeroBias_phi,
               [nmin1,nmax1],[nmin2,nmax2],np.linspace(0,1,1000),dr_min,dr_max)

# %%



