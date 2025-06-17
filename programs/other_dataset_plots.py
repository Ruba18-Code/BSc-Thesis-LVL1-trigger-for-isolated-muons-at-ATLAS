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

nmin=0
nmax=2000
dr_min=0.2
dr_max=0.6

ZeroBias_eta=energy_cut(MuonTree_ZeroBias["muon_pt"].array(), MuonTree_ZeroBias["muon_eta"].array())[nmin:nmax]
ZeroBias_phi=energy_cut(MuonTree_ZeroBias["muon_pt"].array(), MuonTree_ZeroBias["muon_phi"].array())[nmin:nmax]
ZeroBias_pt=energy_cut(MuonTree_ZeroBias["muon_pt"].array(), MuonTree_ZeroBias["muon_pt"].array())[nmin:nmax]
ZeroBias_e=energy_cut(MuonTree_ZeroBias["muon_pt"].array(), MuonTree_ZeroBias["muon_e"].array())[nmin:nmax]
ZeroBias_charge=energy_cut(MuonTree_ZeroBias["muon_pt"].array(), MuonTree_ZeroBias["muon_charge"].array())[nmin:nmax]

#Select quality 0 Z->mumu
Zmumu_pt=quality_selector(MuonTree_Zmumu["muon_quality"].array(),MuonTree_Zmumu["muon_pt"].array(),0)[nmin:nmax]
Zmumu_eta=quality_selector(MuonTree_Zmumu["muon_quality"].array(),MuonTree_Zmumu["muon_eta"].array(),0)[nmin:nmax]
Zmumu_phi=quality_selector(MuonTree_Zmumu["muon_quality"].array(),MuonTree_Zmumu["muon_phi"].array(),0)[nmin:nmax]

#And select the Z peak pairs
Zmumu_pt, Zmumu_eta, Zmumu_phi = get_all_Z_peak_pairs(Zmumu_pt,Zmumu_eta,Zmumu_phi)


# %%
coolplot([ZeroBias_pt, Zmumu_pt], np.linspace(0,2*10**5, 50))

# %%
res1=muon_isolation_all_events(MuonTree_ZeroBias, ZeroBias_eta, ZeroBias_phi, dr_min, dr_max, [nmin,nmax])
res2=muon_isolation_all_events(MuonTree_Zmumu, Zmumu_eta, Zmumu_phi, dr_min, dr_max, [nmin,nmax])

# %%
coolplot([res1,res2], np.linspace(0,10000,50))

# %%
ratio1=ak.flatten(res1)/ak.flatten(ZeroBias_pt)
ratio2=ak.flatten(res2)/ak.flatten(Zmumu_pt)

coolplot([ratio1,ratio2], np.linspace(0,0.4,50))

# %%
dr_min=[0.2,0.3,0.1,0.04]
dr_max=[0.8,0.6,0.4,0.26]

plot_ROC_curve(MuonTree_Zmumu,MuonTree_ZeroBias,Zmumu_pt,Zmumu_eta,Zmumu_phi,ZeroBias_pt,ZeroBias_eta,ZeroBias_phi,
               [nmin,nmax],[nmin,nmax],np.linspace(0,1,5*int(np.sqrt(nmax-nmin))),dr_min,dr_max)

# %%



