import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from my_functions import*


file= uproot.open("/home/ruben/Escritorio/BachelorThesisRuben/Data/Muon_trees.root") #opening the Root file with Uproot 

MuonTree_Zmumu=file["MuonTree_Zmumu;1"]

Zmumu_data=MuonTree_Zmumu["muon_e"].array()

MuonTree_Zmumu=file["MuonTree_Zmumu;1"] 

eta=MuonTree_Zmumu["muon_eta"].array()
phi=MuonTree_Zmumu["muon_phi"].array()
pt=MuonTree_Zmumu["muon_pt"].array()

Z_mass=invariant_mass_all_muons(pt,eta,phi) #computes the invariant mass vector


data=[Z_mass]
bins=np.linspace(6*10**4,1.2*10**5)
colors=['r','b']
x_label=r"Energy (MeV)"
y_label="Counts"
title=r"Z boson invariant mass reconstruction"
label2='Z boson invariant mass'
labels=[label2]
coolplot(data,bins,colors,labels,x_label,y_label,title)


#Let's also plot for high quality data
data1=MuonTree_Zmumu["muon_quality"].array() #save the quality vector
value=0 #set desired quality to 0
data2=MuonTree_Zmumu["muon_pt"].array()
pt=quality_selector(data1,data2,value)
data2=MuonTree_Zmumu["muon_eta"].array()
eta=quality_selector(data1,data2,value)
data2=MuonTree_Zmumu["muon_phi"].array()
phi=quality_selector(data1,data2,value)

Z_mass=invariant_mass_all_muons(pt,eta,phi) #computes the invariant mass vector


data=[Z_mass]
bins=np.linspace(6*10**4,1.2*10**5)
colors=['r','b']
x_label=r"Energy (MeV)"
y_label="Counts"
title=r"Z boson invariant mass reconstruction (only quality 0 data)"
label2='Z boson invariant mass'
labels=[label2]
coolplot(data,bins,colors,labels,x_label,y_label,title)