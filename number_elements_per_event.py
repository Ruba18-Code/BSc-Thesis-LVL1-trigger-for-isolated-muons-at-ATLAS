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

# %%
file= uproot.open("/home/ruben/Escritorio/BachelorThesisRuben/Data/Muon_trees.root") #opening the Root file with Uproot 
MuonTree_Zmumu=file["MuonTree_Zmumu;1"]
MuonTree_ZeroBias=file["MuonTree_ZeroBias;1"]
muon_pt=MuonTree_Zmumu["muon_pt"].array()
muon_pt_back=MuonTree_ZeroBias["muon_e"].array()


# %%
data1=MuonTree_Zmumu["muon_pt"].array()
data1_quality=MuonTree_Zmumu["muon_quality"].array()

data2=MuonTree_ZeroBias["muon_pt"].array()
data2_quality=MuonTree_ZeroBias["muon_quality"].array()


x_label=r"Muon transverse momentum $p_T$ (MeV)"
y_label="Counts"
title=r"Muon counts vs $p_T$ ATLAS Detector (Quality=0)"
label1=r'Z $\longrightarrow \mu \mu$ data'
label2='Zero Bias data'

labels=[label1,label2]

# %%
data=[data1]
bins=np.linspace(0,15,75)
colors=['r','b']
labels=[r'Z $\longrightarrow \mu \mu$ data']
x_label="Number of muons per event"
y_label="Counts"

title="Histogram: number of muons per event"
plot_number_elements_per_event(data,bins,colors,labels,x_label,y_label,title)

# %%
data=[data1,data2]
bins=np.linspace(0,15,75)
colors=['r','b']
labels=[r'Z $\longrightarrow \mu \mu$ data',"Zero Bias data"]
x_label="Number of muons per event"
y_label="Counts"



title="Histogram: number of muons per event"
plot_number_elements_per_event(data,bins,colors,labels,x_label,y_label,title)

# %%
value=0
data1=MuonTree_Zmumu["muon_pt"].array()
data1_quality=MuonTree_Zmumu["muon_quality"].array()
data1=quality_selector(data1_quality,data1,value)

data2=MuonTree_ZeroBias["muon_pt"].array()


data=[data1,data2]
bins=np.linspace(0,15,75)
colors=['r','b']
labels=[r'Z $\longrightarrow \mu \mu$ data',"Zero Bias data"]
x_label="Number of muons per event"
y_label="Counts"
title="Histogram: number of muons per event (only quality 0 events)"

plot_number_elements_per_event(data,bins,colors,labels,x_label,y_label,title)

# %%



