# %%
from my_functions import*

# %%
file= uproot.open("/home/ruben/Escritorio/BachelorThesisRuben/Data/Muon_trees.root") #opening the Root file with Uproot 
MuonTree_ZeroBias=file["MuonTree_ZeroBias;1"]
MuonTree_Zmumu=file["MuonTree_Zmumu;1"]

# %%
nmin1=0
nmax1=100

nmin2=0
nmax2=3000
dr_min=0.0
dr_max=0.4


# %%
dr_min = np.linspace(0, 0, 3)
dr_max = np.linspace(0.2, 2, 3)
bins = np.linspace(0, 5, 120)

plot_ROC_curve(MuonTree_Zmumu, MuonTree_ZeroBias, [0,50],[1000,2000],bins,dr_min,dr_max)

# %%
dr_min = np.linspace(0, 0, 3)
dr_max = np.linspace(0.2, 1, 3)
bins = np.linspace(0, 5, 120)

plot_ROC_curve(MuonTree_Zmumu, MuonTree_ZeroBias, [0,100],[1000,3000],bins,dr_min,dr_max)

# %%
dr_min = np.linspace(0, 0, 3)
dr_max = np.linspace(0.2, 0.5, 3)
bins = np.linspace(0, 5, 120)

plot_ROC_curve(MuonTree_Zmumu, MuonTree_ZeroBias, [0,100],[1000,3000],bins,dr_min,dr_max)

# %%
dr_min = np.linspace(0, 0.3, 3)
dr_max = np.linspace(0.5, 0.5, 3)
bins = np.linspace(0, 5, 120)

plot_ROC_curve(MuonTree_Zmumu, MuonTree_ZeroBias, [0,100],[1000,3000],bins,dr_min,dr_max)

# %%



