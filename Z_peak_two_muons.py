# %%
import my_functions
import uproot
import numpy as np #math and science package
import scipy as sp #math and science package
import awkward as ak #root files are usuallt awkward arrays 
import matplotlib.pyplot as plt #plot stuff
from my_functions import coolplot
from my_functions import ak_element_lenght_counter
from my_functions import plot_number_elements_per_event
from my_functions import quality_locator
from my_functions import quality_selector
from my_functions import ak_element_lenght_counter
from my_functions import invariant_mass_all_muons
from my_functions import invariant_mass

# %%
file= uproot.open("/home/ruben/Escritorio/BachelorThesisRuben/Data/Muon_trees.root") #opening the Root file with Uproot 

MuonTree_Zmumu=file["MuonTree_Zmumu;1"]

Zmumu_data=MuonTree_Zmumu["muon_e"].array()

MuonTree_Zmumu=file["MuonTree_Zmumu;1"] 

eta=MuonTree_Zmumu["muon_eta"].array()
phi=MuonTree_Zmumu["muon_phi"].array()
pt=MuonTree_Zmumu["muon_pt"].array()

Z_mass=invariant_mass_all_muons(pt,eta,phi) #computes the invariant mass vector


# %%
data=[Z_mass]
bins=np.linspace(6*10**4,1.2*10**5)
colors=['r','b']
x_label=r"Energy (MeV)"
y_label="Counts"
title=r"Z boson invariant mass reconstruction"
label2='Z boson invariant mass'
labels=[label2]
coolplot(data,bins,colors,labels,x_label,y_label,title)

# %%
#Now let's make the plot using quality=0 data.

file= uproot.open("/home/ruben/Escritorio/BachelorThesisRuben/Data/Muon_trees.root") #opening the Root file with Uproot 
MuonTree_Zmumu=file["MuonTree_Zmumu;1"]

data1=MuonTree_Zmumu["muon_quality"].array() #save the quality vector
value=0 #set desired quality to 0

#Select quality 0 data for pt, eta and phi

data2=MuonTree_Zmumu["muon_pt"].array()
pt=quality_selector(data1,data2,value)
data2=MuonTree_Zmumu["muon_eta"].array()
eta=quality_selector(data1,data2,value)
data2=MuonTree_Zmumu["muon_phi"].array()
phi=quality_selector(data1,data2,value)

#Compute the Z peak for the quality 0 data

Z_mass_high_quality=invariant_mass_all_muons(pt,eta,phi)

#Plot the Z peak for quality 0 compared with the Z peak using all qualities

data=[Z_mass_high_quality]
bins=np.linspace(6*10**4,1.2*10**5)
colors=['r','b']
x_label=r"Energy (MeV)"
y_label="Counts"
title=r"Z boson invariant mass reconstruction (only quality 0 events)"
label2='Z boson invariant mass'
labels=[label2]

coolplot(data,bins,colors,labels,x_label,y_label,title)

# %%
data=[Z_mass_high_quality, Z_mass]
bins=np.linspace(6*10**4,1.2*10**5)
colors=['r','b']
x_label=r"Energy (MeV)"
y_label="Counts"
title=r"Z boson invariant mass reconstruction"
label1='0 quality'
label2='All quality'
labels=[label1,label2]

coolplot(data,bins,colors,labels,x_label,y_label,title)

# %%
