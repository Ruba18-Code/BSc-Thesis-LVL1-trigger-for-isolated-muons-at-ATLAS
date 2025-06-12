# %%
from my_functions import*

# %%
file= uproot.open("/home/ruben/Escritorio/BachelorThesisRuben/Data/Muon_trees.root") #opening the Root file with Uproot 
MuonTree_ZeroBias=file["MuonTree_ZeroBias;1"]
MuonTree_Zmumu=file["MuonTree_Zmumu;1"]

# %%
nmin1=0
nmax1=1000

nmin2=0
nmax2=1000

Zmumu_pt=quality_selector(MuonTree_Zmumu["muon_quality"].array(),MuonTree_Zmumu["muon_pt"].array(),0)[nmin1:nmax1]
Zmumu_eta=quality_selector(MuonTree_Zmumu["muon_quality"].array(),MuonTree_Zmumu["muon_eta"].array(),0)[nmin1:nmax1]
Zmumu_phi=quality_selector(MuonTree_Zmumu["muon_quality"].array(),MuonTree_Zmumu["muon_phi"].array(),0)[nmin1:nmax1]

ZeroBias_pt=energy_cut(MuonTree_ZeroBias["muon_pt"].array(), MuonTree_ZeroBias["muon_pt"].array(), 14*10**3)[nmin2:nmax2]
ZeroBias_eta=energy_cut(MuonTree_ZeroBias["muon_pt"].array(), MuonTree_ZeroBias["muon_eta"].array(), 14*10**3)[nmin2:nmax2]
ZeroBias_phi=energy_cut(MuonTree_ZeroBias["muon_pt"].array(), MuonTree_ZeroBias["muon_phi"].array(), 14*10**3)[nmin2:nmax2]
# %%
dr_min = np.linspace(0, 0, 3)
dr_max = np.linspace(0.2, 2, 3)
bins = np.linspace(0, 2, 200)

plot_ROC_curve(MuonTree_Zmumu, MuonTree_ZeroBias, Zmumu_pt, Zmumu_eta, Zmumu_phi, ZeroBias_pt, ZeroBias_eta,
                ZeroBias_phi, [nmin1,nmax1],[nmin2,nmax2],bins,dr_min,dr_max,)

# %%
dr_min = np.linspace(0, 0, 3)
dr_max = np.linspace(0.2, 1, 3)
bins = np.linspace(0, 2, 200)

plot_ROC_curve(MuonTree_Zmumu, MuonTree_ZeroBias, Zmumu_pt, Zmumu_eta, Zmumu_phi, ZeroBias_pt, ZeroBias_eta,
                ZeroBias_phi, [nmin1,nmax1],[nmin2,nmax2],bins,dr_min,dr_max,)

# %%
dr_min = np.linspace(0, 0, 3)
dr_max = np.linspace(0.2, 0.5, 3)
bins = np.linspace(0, 2, 200)

plot_ROC_curve(MuonTree_Zmumu, MuonTree_ZeroBias, Zmumu_pt, Zmumu_eta, Zmumu_phi, ZeroBias_pt, ZeroBias_eta,
                ZeroBias_phi, [nmin1,nmax1],[nmin2,nmax2],bins,dr_min,dr_max,)

# %%
dr_min = np.linspace(0, 0.3, 3)
dr_max = np.linspace(0.5, 0.5, 3)
bins = np.linspace(0, 2, 200)

plot_ROC_curve(MuonTree_Zmumu, MuonTree_ZeroBias, Zmumu_pt, Zmumu_eta, Zmumu_phi, ZeroBias_pt, ZeroBias_eta,
                ZeroBias_phi, [nmin1,nmax1],[nmin2,nmax2],bins,dr_min,dr_max,)

# %%



