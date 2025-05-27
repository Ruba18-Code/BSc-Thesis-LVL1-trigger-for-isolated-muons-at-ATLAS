# %%
from my_functions import *

# %%
file= uproot.open("/home/ruben/Escritorio/BachelorThesisRuben/Data/Muon_trees.root") #opening the Root file with Uproot 
MuonTree_ZeroBias=file["MuonTree_ZeroBias;1"]
MuonTree_Zmumu=file["MuonTree_Zmumu;1"]

# %%
muon_eta_all=MuonTree_ZeroBias["LVL1Muon_eta"].array()
muon_phi_all=MuonTree_ZeroBias["LVL1Muon_phi"].array()

non_empty_count = ak.sum(ak.num(muon_eta_all) > 0)
print("Non-empty elements:", non_empty_count)

res1=muon_isolation_all_events(MuonTree_ZeroBias,muon_eta_all,muon_phi_all,0.4,[0,10000],5000)
# %%
muon_eta_all=MuonTree_Zmumu["LVL1Muon_eta"].array()
muon_phi_all=MuonTree_Zmumu["LVL1Muon_phi"].array()

non_empty_count = ak.sum(ak.num(muon_eta_all) > 0)
print("Non-empty elements:", non_empty_count)

res2=res=muon_isolation_all_events(MuonTree_Zmumu,muon_eta_all,muon_phi_all,0.4,[0,1000],500)

# %%
MuonTree_Zmumu=file["MuonTree_Zmumu;1"]

data1=ak.flatten(res1)
data2=ak.flatten(res2)

colors=['b','r']
labels=[r"Zero Bias LVL1 muons",r"Z $\longrightarrow \mu \mu$ LVL1 muons"]

coolplot([data1,data2],np.linspace(0,30000,25),colors,labels,"Energy (MeV)","Counts","Energy histogram isolated muons ATLAS detector")

print(r"Some of the isolated energies for Z to mumu LVL1 muons are:", res2[0] ,"and", res[1], ". The mean "
"value is", ak.mean(res2), "MeV")