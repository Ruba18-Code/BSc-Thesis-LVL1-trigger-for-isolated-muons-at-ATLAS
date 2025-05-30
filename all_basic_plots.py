# %%
import uproot
import numpy as np #math and science package
import scipy as sp #math and science package
import awkward as ak #root files are usuallt awkward arrays 
import matplotlib.pyplot as plt #plot stuff

from my_functions import*
# %%
file= uproot.open("/home/ruben/Escritorio/BachelorThesisRuben/Data/Muon_trees.root") #opening the Root file with Uproot 

file.keys() #Here we can see the keys of the file (index)

# %%
MuonTree_Zmumu=file["MuonTree_Zmumu;1"] 
MuonTree_Back=file["MuonTree_ZeroBias;1"]

MuonTree_Zmumu.show() #Let's see which plots can we do. We're going to focus on the offline muons (the ones with name muon_something).
                      #We're going to plot the Zmumu data together with the 0 bias (background) in order to compare them.

# %%
data1=MuonTree_Zmumu["muon_pt"].array()
data2=MuonTree_Back["muon_pt"].array()

data=[data1,data2]
bins=np.linspace(0,10**5,50)
colors=['r','b']
x_label=r"Muon transverse momentum $p_T$ (MeV)"
y_label="Counts"
title=r"Muon counts vs $p_T$ ATLAS Detector"
label1=r'Z $\longrightarrow \mu \mu$ data'
label2='Background data'
labels=[label1,label2]

coolplot(data,bins,colors,labels,x_label,y_label,title)


# %%
data1=MuonTree_Zmumu["muon_eta"].array()
data2=MuonTree_Back["muon_eta"].array()

data=[data1,data2]
bins=np.linspace(-4,4,50)
colors=['r','b']
x_label=r"Muon pseudorapitidy $\eta$"
y_label="Counts"
title=r"Muon counts vs $\eta$ ATLAS Detector"
label1=r'Z $\longrightarrow \mu \mu$ data'
label2='Background data'
labels=[label1,label2]

coolplot(data,bins,colors,labels,x_label,y_label,title)

# %%
data1=MuonTree_Zmumu["muon_phi"].array()
data2=MuonTree_Back["muon_phi"].array()

data=[data1,data2]
bins=np.linspace(-4,4,50)
colors=['r','b']
x_label=r"Muon phi $\phi$"
y_label="Counts"
title=r"Muon counts vs $\phi$ ATLAS Detector"
label1=r'Z $\longrightarrow \mu \mu$ data'
label2='Background data'
labels=[label1,label2]

coolplot(data,bins,colors,labels,x_label,y_label,title)

# %%
data1=MuonTree_Zmumu["muon_e"].array()
data2=MuonTree_Back["muon_e"].array()

data=[data1,data2]
bins=np.linspace(0,2*10**5,50)
colors=['r','b']
x_label=r"Muon Energy $E$ (MeV)"
y_label="Counts"
title=r"Muon counts vs $E$ ATLAS Detector"
label1=r'Z $\longrightarrow \mu \mu$ data'
label2='Background data'
labels=[label1,label2]

coolplot(data,bins,colors,labels,x_label,y_label,title)

# %%
data1=MuonTree_Zmumu["muon_type"].array()
data2=MuonTree_Back["muon_type"].array()

data=[data1,data2]
bins=np.linspace(5,10,50)
colors=['r','b']
x_label=r"Muon type"
y_label="Counts"
title=r"Muon counts vs type ATLAS Detector"
label1=r'Z $\longrightarrow \mu \mu$ data'
label2='Background data'
labels=[label1,label2]

coolplot(data,bins,colors,labels,x_label,y_label,title)

# %%
data1=MuonTree_Zmumu["muon_charge"].array()
data2=MuonTree_Back["muon_charge"].array()

data=[data1,data2]
bins=np.linspace(-2,2,50)
colors=['r','b']
x_label=r"Muon charge (e-)"
y_label="Counts"
title=r"Muon counts vs charge ATLAS Detector"
label1=r'Z $\longrightarrow \mu \mu$ data'
label2='Background data'
labels=[label1,label2]

coolplot(data,bins,colors,labels,x_label,y_label,title)

# %%
data1=MuonTree_Zmumu["muon_author"].array()
data2=MuonTree_Back["muon_author"].array()

data=[data1,data2]
bins=np.linspace(0,12,50)
colors=['r','b']
x_label=r"Muon author"
y_label="Counts"
title=r"Muon counts vs author ATLAS Detector"
label1=r'Z $\longrightarrow \mu \mu$ data'
label2='Background data'
labels=[label1,label2]

coolplot(data,bins,colors,labels,x_label,y_label,title)

# %%
data1=MuonTree_Zmumu["muon_quality"].array()
data2=MuonTree_Back["muon_quality"].array()

data=[data1,data2]
bins=np.linspace(-1,10,50)
colors=['r','b']
x_label=r"Muon quality"
y_label="Counts"
title=r"Muon counts vs quality ATLAS Detector"
label1=r'Z $\longrightarrow \mu \mu$ data'
label2='Background data'
labels=[label1,label2]

coolplot(data,bins,colors,labels,x_label,y_label,title)

# %% [markdown]
# ## QUALITY 0 PLOTS

# %%
value=0 

MuonTree_Zmumu=file["MuonTree_Zmumu;1"] 
MuonTree_ZeroBias=file["MuonTree_ZeroBias;1"]

# %%
#Select the data and keep only 0 quality values

data1=MuonTree_Zmumu["muon_pt"].array()
data1_quality=MuonTree_Zmumu["muon_quality"].array()
data1=quality_selector(data1_quality,data1,0)

#For the Zero bias data we do not select quality (because it's noise anyways)

data2=MuonTree_ZeroBias["muon_pt"].array()
data2_quality=MuonTree_ZeroBias["muon_quality"].array()

data=[data1,data2]
bins=np.linspace(0,10**5,50)
colors=['r','b']
x_label=r"Muon transverse momentum $p_T$ (MeV)"
y_label="Counts"
title=r"Muon counts vs $p_T$ ATLAS Detector (quality 0 data)"
label1=r'Z $\longrightarrow \mu \mu$ data'
label2='Background data'
labels=[label1,label2]

coolplot(data,bins,colors,labels,x_label,y_label,title)

# %%
#Select the data and keep only 0 quality values

data1=MuonTree_Zmumu["muon_eta"].array()
data1_quality=MuonTree_Zmumu["muon_quality"].array()
data1=quality_selector(data1_quality,data1,0)

#For the Zero bias data we do not select quality (because it's noise anyways)

data2=MuonTree_ZeroBias["muon_eta"].array()
data2_quality=MuonTree_ZeroBias["muon_quality"].array()

data=[data1,data2]
bins=np.linspace(-4,4,50)
colors=['r','b']
x_label=r"Muon pseudorapidity $\eta$"
y_label="Counts"
title=r"Muon counts vs $\eta$ ATLAS Detector (quality 0 data)"
label1=r'Z $\longrightarrow \mu \mu$ data'
label2='Background data'
labels=[label1,label2]

coolplot(data,bins,colors,labels,x_label,y_label,title)

# %%
#Select the data and keep only 0 quality values

data1=MuonTree_Zmumu["muon_phi"].array()
data1_quality=MuonTree_Zmumu["muon_quality"].array()
data1=quality_selector(data1_quality,data1,0)

#For the Zero bias data we do not select quality (because it's noise anyways)

data2=MuonTree_ZeroBias["muon_phi"].array()
data2_quality=MuonTree_ZeroBias["muon_quality"].array()

data=[data1,data2]
bins=np.linspace(-4,4,50)
colors=['r','b']
x_label=r"Muon phi $\phi$"
y_label="Counts"
title=r"Muon counts vs $\phi$ ATLAS Detector (quality 0 data)"
label1=r'Z $\longrightarrow \mu \mu$ data'
label2='Background data'
labels=[label1,label2]

coolplot(data,bins,colors,labels,x_label,y_label,title)

# %%
#Select the data and keep only 0 quality values

data1=MuonTree_Zmumu["muon_e"].array()
data1_quality=MuonTree_Zmumu["muon_quality"].array()
data1=quality_selector(data1_quality,data1,0)

#For the Zero bias data we do not select quality (because it's noise anyways)

data2=MuonTree_ZeroBias["muon_e"].array()
data2_quality=MuonTree_ZeroBias["muon_quality"].array()

data=[data1,data2]
bins=np.linspace(0,2*10**5,50)
colors=['r','b']
x_label=r"Muon energy $E$ (MeV)"
y_label="Counts"
title=r"Muon counts vs $E$ ATLAS Detector (quality 0 data)"
label1=r'Z $\longrightarrow \mu \mu$ data'
label2='Background data'
labels=[label1,label2]

coolplot(data,bins,colors,labels,x_label,y_label,title)

# %%
#Select the data and keep only 0 quality values

data1=MuonTree_Zmumu["muon_type"].array()
data1_quality=MuonTree_Zmumu["muon_quality"].array()
data1=quality_selector(data1_quality,data1,0)

#For the Zero bias data we do not select quality (because it's noise anyways)

data2=MuonTree_ZeroBias["muon_type"].array()
data2_quality=MuonTree_ZeroBias["muon_quality"].array()

data=[data1,data2]
bins=np.linspace(0,10,50)
colors=['r','b']
x_label=r"Muon type"
y_label="Counts"
title=r"Muon counts vs type ATLAS Detector (quality 0 data)"
label1=r'Z $\longrightarrow \mu \mu$ data'
label2='Background data'
labels=[label1,label2]

coolplot(data,bins,colors,labels,x_label,y_label,title)

# %%
#Select the data and keep only 0 quality values

data1=MuonTree_Zmumu["muon_charge"].array()
data1_quality=MuonTree_Zmumu["muon_quality"].array()
data1=quality_selector(data1_quality,data1,0)

#For the Zero bias data we do not select quality (because it's noise anyways)

data2=MuonTree_ZeroBias["muon_charge"].array()
data2_quality=MuonTree_ZeroBias["muon_quality"].array()

data=[data1,data2]
bins=np.linspace(-2,2,50)
colors=['r','b']
x_label=r"Muon charge (e-)"
y_label="Counts"
title=r"Muon counts vs charge ATLAS Detector (quality 0 data)"
label1=r'Z $\longrightarrow \mu \mu$ data'
label2='Background data'
labels=[label1,label2]

coolplot(data,bins,colors,labels,x_label,y_label,title)

# %%
#Select the data and keep only 0 quality values

data1=MuonTree_Zmumu["muon_author"].array()
data1_quality=MuonTree_Zmumu["muon_quality"].array()
data1=quality_selector(data1_quality,data1,0)

#For the Zero bias data we do not select quality (because it's noise anyways)

data2=MuonTree_ZeroBias["muon_author"].array()
data2_quality=MuonTree_ZeroBias["muon_quality"].array()

data=[data1,data2]
bins=np.linspace(0,10,50)
colors=['r','b']
x_label=r"Muon author"
y_label="Counts"
title=r"Muon counts vs author ATLAS Detector (quality 0 data)"
label1=r'Z $\longrightarrow \mu \mu$ data'
label2='Background data'
labels=[label1,label2]

coolplot(data,bins,colors,labels,x_label,y_label,title)

# %%
#Select the data and keep only 0 quality values

data1=MuonTree_Zmumu["muon_quality"].array()
data1_quality=MuonTree_Zmumu["muon_quality"].array()
data1=quality_selector(data1_quality,data1,0)

#For the Zero bias data we do not select quality (because it's noise anyways)

data2=MuonTree_ZeroBias["muon_quality"].array()
data2_quality=MuonTree_ZeroBias["muon_quality"].array()

data=[data1,data2]
bins=np.linspace(-1,10,50)
colors=['r','b']
x_label=r"Muon quality"
y_label="Counts"
title=r"Muon counts vs quality ATLAS Detector (quality 0 data)"
label1=r'Z $\longrightarrow \mu \mu$ data'
label2='Background data'
labels=[label1,label2]

coolplot(data,bins,colors,labels,x_label,y_label,title)


