# %%
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from my_functions import*

# %%
file= uproot.open("/home/ruben/Escritorio/BachelorThesisRuben/Data/Muon_trees.root") #opening the Root file with Uproot 
MuonTree_ZeroBias=file["MuonTree_ZeroBias;1"]
MuonTree_Zmumu=file["MuonTree_Zmumu;1"]

# %%
nmin1=5000
nmax1=7000

nmin2=5000
nmax2=7000

Zmumu_pt=quality_selector(MuonTree_Zmumu["muon_quality"].array(),MuonTree_Zmumu["muon_pt"].array(),0)[nmin1:nmax1]
Zmumu_eta=quality_selector(MuonTree_Zmumu["muon_quality"].array(),MuonTree_Zmumu["muon_eta"].array(),0)[nmin1:nmax1]
Zmumu_phi=quality_selector(MuonTree_Zmumu["muon_quality"].array(),MuonTree_Zmumu["muon_phi"].array(),0)[nmin1:nmax1]

ZeroBias_pt=energy_cut(MuonTree_ZeroBias["muon_pt"].array(), MuonTree_ZeroBias["muon_pt"].array(), 14*10**3)[nmin2:nmax2]
ZeroBias_eta=energy_cut(MuonTree_ZeroBias["muon_pt"].array(), MuonTree_ZeroBias["muon_eta"].array(), 14*10**3)[nmin2:nmax2]
ZeroBias_phi=energy_cut(MuonTree_ZeroBias["muon_pt"].array(), MuonTree_ZeroBias["muon_phi"].array(), 14*10**3)[nmin2:nmax2]
# %%
dr_min=[0,0,0,0]
dr_max=[0.25,0.4,0.6,0.8]
bins = np.linspace(0, 1, 1000)

plot_ROC_curve(MuonTree_Zmumu, MuonTree_ZeroBias, Zmumu_pt, Zmumu_eta, Zmumu_phi, ZeroBias_pt, ZeroBias_eta,
                ZeroBias_phi, [nmin1,nmax1],[nmin2,nmax2],bins,dr_min,dr_max,)

# %%
dr_min = np.linspace(0, 0, 3)
dr_max = np.linspace(0.2, 1, 3)
bins = np.linspace(0, 1, 1000)

plot_ROC_curve(MuonTree_Zmumu, MuonTree_ZeroBias, Zmumu_pt, Zmumu_eta, Zmumu_phi, ZeroBias_pt, ZeroBias_eta,
                ZeroBias_phi, [nmin1,nmax1],[nmin2,nmax2],bins,dr_min,dr_max,)

# %%
dr_min = np.linspace(0, 0, 3)
dr_max = np.linspace(0.2, 0.5, 3)
bins = np.linspace(0, 1, 1000)

plot_ROC_curve(MuonTree_Zmumu, MuonTree_ZeroBias, Zmumu_pt, Zmumu_eta, Zmumu_phi, ZeroBias_pt, ZeroBias_eta,
                ZeroBias_phi, [nmin1,nmax1],[nmin2,nmax2],bins,dr_min,dr_max,)

dr_min=[0,0.05,0.1,0.05]
dr_max=[0.3,0.3,0.4,0.25]
bins = np.linspace(0, 2, 2000)

plot_ROC_curve(MuonTree_Zmumu, MuonTree_ZeroBias, Zmumu_pt, Zmumu_eta, Zmumu_phi, ZeroBias_pt, ZeroBias_eta,
                ZeroBias_phi, [nmin1,nmax1],[nmin2,nmax2],bins,dr_min,dr_max,)

# %%



