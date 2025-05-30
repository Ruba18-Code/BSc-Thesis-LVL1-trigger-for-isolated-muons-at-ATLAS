# %%
# %%
from my_functions import *

# %%
file= uproot.open("/home/ruben/Escritorio/BachelorThesisRuben/Data/Muon_trees.root") #opening the Root file with Uproot 
MuonTree_ZeroBias=file["MuonTree_ZeroBias;1"]
MuonTree_Zmumu=file["MuonTree_Zmumu;1"]

# %%
"""
Here I'm going to use LVL1 muons
"""

muon_eta_all=MuonTree_ZeroBias["LVL1Muon_eta"].array()
muon_phi_all=MuonTree_ZeroBias["LVL1Muon_phi"].array()

non_empty_count = ak.sum(ak.num(muon_eta_all) > 0)
print("Non-empty elements:", non_empty_count)

res1=muon_isolation_all_events(MuonTree_ZeroBias,muon_eta_all,muon_phi_all,0.05,0.4,[0,5000],500)
# %%
muon_eta_all=MuonTree_Zmumu["LVL1Muon_eta"].array()
muon_phi_all=MuonTree_Zmumu["LVL1Muon_phi"].array()

non_empty_count = ak.sum(ak.num(muon_eta_all) > 0)
print("Non-empty elements:", non_empty_count)

res2=res=muon_isolation_all_events(MuonTree_Zmumu,muon_eta_all,muon_phi_all,0.05,0.4,[0,1000],100)

# %%
data1=ak.flatten(res1)
data2=ak.flatten(res2)

colors=['b','r']
labels=[r"Zero Bias LVL1 muons",r"Z $\longrightarrow \mu \mu$ LVL1 muons"]

coolplot([data1,data2],np.linspace(0,30000,25),colors,labels,"Transverse energy (MeV)","Counts","Energy histogram isolated muons ATLAS detector")

print(r"Some of the isolated energies for Z to mumu LVL1 muons are:", res2[0] ,"and", res[1], ". The mean "
"value is", ak.mean(res2), "MeV")


# %%

"""
Here I want to use offline muons
"""

muon_eta_all=MuonTree_ZeroBias["muon_eta"].array()
muon_phi_all=MuonTree_ZeroBias["muon_phi"].array()

non_empty_count = ak.sum(ak.num(muon_eta_all) > 0)
print("Non-empty elements:", non_empty_count)

res1=muon_isolation_all_events(MuonTree_ZeroBias,muon_eta_all,muon_phi_all,0.05,0.4,[12000,14000],500)
# %%
muon_eta_all=MuonTree_Zmumu["muon_eta"].array()
muon_phi_all=MuonTree_Zmumu["muon_phi"].array()

non_empty_count = ak.sum(ak.num(muon_eta_all) > 0)
print("Non-empty elements:", non_empty_count)

res2=res=muon_isolation_all_events(MuonTree_Zmumu,muon_eta_all,muon_phi_all,0.05,0.4,[0,1000],500)

data1=ak.flatten(res1)
data2=ak.flatten(res2)

colors=['b','r']
labels=[r"Zero Bias offline muons",r"Z $\longrightarrow \mu \mu$ offline muons"]

coolplot([data1,data2],np.linspace(0,30000,25),colors,labels,"Transverse energy (MeV)","Counts","Energy histogram isolated muons ATLAS detector")

print(r"Some of the isolated energies for Z to mumu offline muons are:", res2[0] ,"and", res[1], ". The mean "
"value is", ak.mean(res2), "MeV")

############################################################


# %%
"""
Here I want to plot the same thing but I'm going to select only quality 0 data for the Z->mu mu events
"""

muon_eta_all=MuonTree_ZeroBias["muon_eta"].array()
muon_phi_all=MuonTree_ZeroBias["muon_phi"].array()

non_empty_count = ak.sum(ak.num(muon_eta_all) > 0)
print("Non-empty elements:", non_empty_count)

res1=muon_isolation_all_events(MuonTree_ZeroBias,muon_eta_all,muon_phi_all,0.05,0.4,[12000,14000],500)

"""
Now select only quality 0 events for the Z mumu data
"""

muon_eta_all=quality_selector(MuonTree_Zmumu["muon_quality"].array(),MuonTree_Zmumu["muon_eta"].array(),0)
muon_phi_all=quality_selector(MuonTree_Zmumu["muon_quality"].array(),MuonTree_Zmumu["muon_phi"].array(),0)

non_empty_count = ak.sum(ak.num(muon_eta_all) > 0)
print("Non-empty elements:", non_empty_count)

res2=res=muon_isolation_all_events(MuonTree_Zmumu,muon_eta_all,muon_phi_all,0.05,0.4,[0,1000],500)

data1=ak.flatten(res1)
data2=ak.flatten(res2)

colors=['b','r']
labels=[r"Zero Bias offline muons",r"Z $\longrightarrow \mu \mu$ quality 0 offline muons"]

coolplot([data1,data2],np.linspace(0,30000,25),colors,labels,"Transverse energy (MeV)","Counts"
,"Energy histogram isolated muons ATLAS detector- Only quality 0 events")

print(r"Some of the isolated energies for Z to mumu quality 0 offline muons are:", res2[0] ,"and", res2[1], ". The mean "
"value is", ak.mean(res2), "MeV")


