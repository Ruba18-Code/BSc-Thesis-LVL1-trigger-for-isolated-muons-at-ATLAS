import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from my_functions import*
# %%
file= uproot.open("/home/ruben/Escritorio/BachelorThesisRuben/Data/Muon_trees.root") #opening the Root file with Uproot 

file.keys() #Here we can see the keys of the file (index)

# %%
MuonTree_Zmumu=file["MuonTree_Zmumu;1"] 
MuonTree_ZeroBias=file["MuonTree_ZeroBias;1"]

MuonTree_Zmumu.show() #Let's see which plots can we do. We're going to focus on the offline muons (the ones with name muon_something).
                      #We're going to plot the Zmumu data together with the 0 bias (background) in order to compare them.

# %%
#Set the range of events to plot
nmin1=0
nmax1=10000

#Select quality 0 Z->mumu
Zmumu_pt=quality_selector_with_empty(MuonTree_Zmumu["muon_quality"].array(),MuonTree_Zmumu["muon_pt"].array(),0)[nmin1:nmax1]
Zmumu_eta=quality_selector_with_empty(MuonTree_Zmumu["muon_quality"].array(),MuonTree_Zmumu["muon_eta"].array(),0)[nmin1:nmax1]
Zmumu_phi=quality_selector_with_empty(MuonTree_Zmumu["muon_quality"].array(),MuonTree_Zmumu["muon_phi"].array(),0)[nmin1:nmax1]
Zmumu_e=quality_selector_with_empty(MuonTree_Zmumu["muon_quality"].array(),MuonTree_Zmumu["muon_e"].array(),0)[nmin1:nmax1]
Zmumu_charge=quality_selector_with_empty(MuonTree_Zmumu["muon_quality"].array(),MuonTree_Zmumu["muon_charge"].array(),0)[nmin1:nmax1]

#And select the Z peak pairs
Zmumu_pt, Zmumu_eta, Zmumu_phi = get_all_Z_peak_pairs(Zmumu_pt,Zmumu_eta,Zmumu_phi)

nmin2=0
nmax2=140000
#Select the ZeroBias data with energy cut
ZeroBias_eta=energy_cut_with_empty(MuonTree_ZeroBias["muon_pt"].array(), MuonTree_ZeroBias["muon_eta"].array())[nmin2:nmax2]
ZeroBias_phi=energy_cut_with_empty(MuonTree_ZeroBias["muon_pt"].array(), MuonTree_ZeroBias["muon_phi"].array())[nmin2:nmax2]
ZeroBias_pt=energy_cut_with_empty(MuonTree_ZeroBias["muon_pt"].array(), MuonTree_ZeroBias["muon_pt"].array())[nmin2:nmax2]
ZeroBias_e=energy_cut_with_empty(MuonTree_ZeroBias["muon_pt"].array(), MuonTree_ZeroBias["muon_e"].array())[nmin2:nmax2]
ZeroBias_charge=energy_cut_with_empty(MuonTree_ZeroBias["muon_pt"].array(), MuonTree_ZeroBias["muon_charge"].array())[nmin2:nmax2]
data1=Zmumu_pt
data2=ZeroBias_pt
l1=len(ak.flatten(data1))
l2=len(ak.flatten(data2))
data=[data1,data2]
bins=np.linspace(0,1.5*10**5,50)
colors=['r','b']
x_label=r"Muon transverse momentum $p_T$ (MeV)"
y_label="Counts"
title=r"Muon counts vs $p_T$ ATLAS Detector"
label1=rf'Z $\longrightarrow \mu \mu$ data, muons={l1}'
label2=f'Zero Bias data, muons={l2}'
labels=[label1,label2]

coolplot(data,bins,colors,labels,x_label,y_label,title)
# %%
data1=Zmumu_eta
data2=ZeroBias_eta
l1=len(ak.flatten(data1))
l2=len(ak.flatten(data2))

data=[data1,data2]
bins=np.linspace(-4,4,50)
colors=['r','b']
x_label=r"Muon pseudorapitidy $\eta$"
y_label="Counts"
title=r"Muon counts vs $\eta$ ATLAS Detector"
label1=rf'Z $\longrightarrow \mu \mu$ data, muons={l1}'
label2=f'Zero Bias data, muons={l2}'
labels=[label1,label2]

coolplot(data,bins,colors,labels,x_label,y_label,title)

# %%
data1=Zmumu_phi
data2=ZeroBias_phi
l1=len(ak.flatten(data1))
l2=len(ak.flatten(data2))

data=[data1,data2]
bins=np.linspace(-4,4,50)
colors=['r','b']
x_label=r"Muon phi $\phi$"
y_label="Counts"
title=r"Muon counts vs $\phi$ ATLAS Detector"
label1=rf'Z $\longrightarrow \mu \mu$ data, muons={l1}'
label2=f'Zero Bias data, muons={l2}'
labels=[label1,label2]

coolplot(data,bins,colors,labels,x_label,y_label,title)

# %%
data1=Zmumu_e
data2=ZeroBias_e
l1=len(ak.flatten(data1))
l2=len(ak.flatten(data2))

data=[data1,data2]
bins=np.linspace(0,2*10**5,50)
colors=['r','b']
x_label=r"Muon Energy $E$ (MeV)"
y_label="Counts"
title=r"Muon counts vs $E$ ATLAS Detector"
label1=rf'Z $\longrightarrow \mu \mu$ data, muons={l1}'
label2=f'Zero Bias data, muons={l2}'
labels=[label1,label2]

coolplot(data,bins,colors,labels,x_label,y_label,title)

# %%
data1=Zmumu_charge
data2=ZeroBias_charge
l1=len(ak.flatten(data1))
l2=len(ak.flatten(data2))

data=[data1,data2]
bins=np.linspace(-2,2,50)
colors=['r','b']
x_label=r"Muon charge (e-)"
y_label="Counts"
title=r"Muon counts vs charge ATLAS Detector"
label1=rf'Z $\longrightarrow \mu \mu$ data, muons={l1}'
label2=f'Zero Bias data, muons={l2}'
labels=[label1,label2]

coolplot(data,bins,colors,labels,x_label,y_label,title)