# %%
from my_functions import*

# %%
file= uproot.open("/home/ruben/Escritorio/BachelorThesisRuben/Data/Muon_trees.root") #opening the Root file with Uproot 

# %%
MuonTree_ZeroBias=file["MuonTree_ZeroBias;1"]
MuonTree_Zmumu=file["MuonTree_Zmumu;1"]

nmin=0
nmax=2000
#Select quality 0 Z->mumu
Zmumu_pt=quality_selector(MuonTree_Zmumu["muon_quality"].array(),MuonTree_Zmumu["muon_pt"].array(),0)[nmin:nmax]
Zmumu_eta=quality_selector(MuonTree_Zmumu["muon_quality"].array(),MuonTree_Zmumu["muon_eta"].array(),0)[nmin:nmax]
Zmumu_phi=quality_selector(MuonTree_Zmumu["muon_quality"].array(),MuonTree_Zmumu["muon_phi"].array(),0)[nmin:nmax]
#And select the Z peak pairs
Zmumu_pt, Zmumu_eta, Zmumu_phi = get_all_Z_peak_pairs(Zmumu_pt,Zmumu_eta,Zmumu_phi)
#Select the ZeroBias data with energy cut
ZeroBias_eta=energy_cut(MuonTree_ZeroBias["muon_pt"].array(), MuonTree_ZeroBias["muon_eta"].array())[nmin:nmax]
ZeroBias_phi=energy_cut(MuonTree_ZeroBias["muon_pt"].array(), MuonTree_ZeroBias["muon_phi"].array())[nmin:nmax]
ZeroBias_pt=energy_cut(MuonTree_ZeroBias["muon_pt"].array(), MuonTree_ZeroBias["muon_pt"].array())[nmin:nmax]
ZeroBias_e=energy_cut(MuonTree_ZeroBias["muon_pt"].array(), MuonTree_ZeroBias["muon_e"].array())[nmin:nmax]

# %%

#Set dr
dr_min=0.05
dr_max=0.3
#Compute isolation
res1=muon_isolation_all_events(MuonTree_ZeroBias, ZeroBias_eta, ZeroBias_phi, dr_min, dr_max, [nmin,nmax])
res2=muon_isolation_all_events(MuonTree_Zmumu, Zmumu_eta, Zmumu_phi, dr_min, dr_max, [nmin,nmax])
#Compute ratio
data1=ak.flatten(res1)/ak.flatten(ZeroBias_pt)
data2=ak.flatten(res2)/ak.flatten(Zmumu_pt)

# %%
d1=data1[data1/data1 == 1]
d2=data2[data2/data2 == 1]


"""
mask=data1/data1 == 1
print(data1[mask]) 
print(ak.flatten(res1)[mask])
print(ak.flatten(ZeroBias_pt)[mask])

# %%
d1=data1[mask]
for indices in ak.where(data1[mask] < 0.0125):
    print(data1[mask][indices]) 
    print(ak.flatten(res1)[mask][indices])
    print(ak.flatten(ZeroBias_pt)[mask][indices])
    print(ak.flatten(ZeroBias_eta)[mask][indices])
    print(ak.flatten(ZeroBias_phi)[mask][indices])
"""
# %%

plot_ROC_curve(MuonTree_Zmumu, MuonTree_ZeroBias, Zmumu_pt, Zmumu_eta, Zmumu_phi, ZeroBias_pt, ZeroBias_eta, ZeroBias_phi,
               [nmin, nmax],[nmin,nmax], np.linspace(0.02,1,1000), [dr_min], [dr_max], title="ROC curve - removing first bins comparison",
               plot_show=False)
plot_ROC_curve(MuonTree_Zmumu, MuonTree_ZeroBias, Zmumu_pt, Zmumu_eta, Zmumu_phi, ZeroBias_pt, ZeroBias_eta, ZeroBias_phi,
               [nmin, nmax],[nmin,nmax], np.linspace(0,1,1000), [dr_min], [dr_max], title="ROC curve - removing first bins comparison",
               plot_show=False)
plt.show()
# %%



