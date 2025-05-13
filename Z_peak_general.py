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
import itertools
from my_functions import pair_selector
from my_functions import invariant_mass_all_muons2


file= uproot.open("/home/ruben/Escritorio/BachelorThesisRuben/Data/Muon_trees.root") #opening the Root file with Uproot 

MuonTree_Zmumu=file["MuonTree_Zmumu;1"]

Zmumu_data=MuonTree_Zmumu["muon_e"].array()

MuonTree_Zmumu=file["MuonTree_Zmumu;1"] 

eta=MuonTree_Zmumu["muon_eta"].array()
phi=MuonTree_Zmumu["muon_phi"].array()
pt=MuonTree_Zmumu["muon_pt"].array()

Z_mass=invariant_mass_all_muons2(pt,eta,phi) #computes the invariant mass vector


data=[Z_mass]
bins=np.linspace(6*10**4,1.2*10**5)
colors=['r','b']
x_label=r"Energy (MeV)"
y_label="Counts"
title=r"Z boson invariant mass reconstruction (more than 2 muons test)"
label2='Z boson invariant mass'
labels=[label2]
coolplot(data,bins,colors,labels,x_label,y_label,title)
