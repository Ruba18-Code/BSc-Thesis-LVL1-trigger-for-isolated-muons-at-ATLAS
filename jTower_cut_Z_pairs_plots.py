# %%
from my_functions import*

# %%
file= uproot.open("/home/ruben/Escritorio/BachelorThesisRuben/Data/Muon_trees.root") #opening the Root file with Uproot 
MuonTree_ZeroBias=file["MuonTree_ZeroBias;1"]
MuonTree_Zmumu=file["MuonTree_Zmumu;1"]

# %%
#Prepare the data for the plots

#Choose the range of events to plot
nmin=0
nmax=1000

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

# %%
#Check how many events are not empty
non_empty_count1= ak.sum(ak.num(ZeroBias_eta[nmin:nmax]) > 0)
non_empty_count2 = ak.sum(ak.num(Zmumu_eta[nmin:nmax]) > 0)

#Compute the isolation and prepare it for plotting ZeroBias
res=muon_isolation_all_events(MuonTree_ZeroBias,ZeroBias_eta,ZeroBias_phi,0.0,0.4,[nmin,nmax],100)
data1=ak.flatten(res)

#Compute the isolation and prepare it for plotting Z mu mu
res=muon_isolation_all_events(MuonTree_Zmumu,Zmumu_eta,Zmumu_phi,0.0,0.4,[nmin,nmax],100)
data2=ak.flatten(res)

colors=["#0072B2", "#FD0000"]
labels=[fr"Zero Bias (cut=14GeV), events={non_empty_count1}",fr"Z $\longrightarrow \mu \mu$ (Q=0, Z peak pairs), events={non_empty_count2}"]

#Plot the data
coolplot([data1,data2],np.linspace(0,30000,40),colors,labels, "Isolation (MeV)","Counts","Isolation histogram")

# %%
#Prepare the data
ratio1=data1/ak.flatten(ZeroBias_pt)
ratio2=data2/ak.flatten(Zmumu_pt)

#Plot the data
colors=["#0072B2", "#FD0000"]
labels=[fr"Zero Bias (cut=14GeV), events={non_empty_count1}",fr"Z $\longrightarrow \mu \mu$ (Q=0, Z peak pairs), events={non_empty_count2}"]

coolplot([ratio1,ratio2],np.linspace(0,1,40),colors,labels, "Isolation/et ratio","Counts","Isolation/et ratio histogram")

# %%
#Now let's do the ROC curve

#Now let's do the ROC curve
#Number of bins is the square root of the number of events 
#The range starts at -1/sqrt(nmax-nmin), this ensures that the data starts at (0,0) even on extreme cases
bins=np.linspace(0,1,int(np.sqrt(nmax-nmin)))
dr_min=[0.0,0.0,0.0]
dr_max=[0.3,0.4,0.5]

#Prepare the data for the plots
plot_ROC_curve(MuonTree_Zmumu, MuonTree_ZeroBias, Zmumu_pt, Zmumu_eta, Zmumu_phi, ZeroBias_pt, ZeroBias_eta,
                ZeroBias_phi, [nmin,nmax],[nmin,nmax],bins,dr_min,dr_max)


