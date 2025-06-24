# %%

import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from my_functions import *
# %%
file= uproot.open("/home/ruben/Escritorio/BachelorThesisRuben/Data/Muon_trees.root") #opening the Root file with Uproot 

# %%
MuonTree_ZeroBias=file["MuonTree_ZeroBias;1"]
MuonTree_Zmumu=file["MuonTree_Zmumu;1"]

nmin=3000
nmax=6000
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
#------------------1st subplot------------------------------------------------------------------------------------------------
#Set scaling factor
scaling=0.5
#Set dr
dr_min=0.05
dr_max=0.32
#Compute isolation
res1=muon_isolation_all_events(MuonTree_ZeroBias, ZeroBias_eta, ZeroBias_phi, dr_min, dr_max, [nmin,nmax], scaling=scaling)
res2=muon_isolation_all_events(MuonTree_Zmumu, Zmumu_eta, Zmumu_phi, dr_min, dr_max, [nmin,nmax], scaling=scaling)
#Compute ratio
data1=ak.flatten(res1)/ak.flatten(ZeroBias_pt)
data2=ak.flatten(res2)/ak.flatten(Zmumu_pt)
#Compute number of non empty events
events1=ak.sum(ak.num(ZeroBias_pt) > 0)
events2=ak.sum(ak.num(Zmumu_pt) > 0)
#Prepare plot
data_1=[data1, data2]
bins=np.linspace(0.0,0.3,40)
colors=["#0072B2", "#FD0000"]
label1=f"Zero Bias, events={events1}"
label2=fr"Z $\to \mu \mu$, events={events2}"
labels=[label1, label2]
x_label=r"$\left( \frac{E_{iso}}{E_{t,muon}} \right)$"
y_label=r"Counts"
title = rf"$\left( \frac{{E_{{iso}}}}{{E_{{t,muon}}}}\right)$ ratio histogram, scaling={scaling}"
#Plot
fig, axis = plt.subplots(3, 2, figsize=(12, 8))
coolplot(data_1, bins, colors, labels, x_label, y_label, title, ax=axis[0, 0])

#------------------2nd subplot------------------------------------------------------------------------------------------------
#Set scaling factor
scaling=1
#Compute isolation
res1=muon_isolation_all_events(MuonTree_ZeroBias, ZeroBias_eta, ZeroBias_phi, dr_min, dr_max, [nmin,nmax], scaling=scaling)
res2=muon_isolation_all_events(MuonTree_Zmumu, Zmumu_eta, Zmumu_phi, dr_min, dr_max, [nmin,nmax], scaling=scaling)
#Compute ratio
data1=ak.flatten(res1)/ak.flatten(ZeroBias_pt)
data2=ak.flatten(res2)/ak.flatten(Zmumu_pt)
#Compute number of non empty events
events1=ak.sum(ak.num(ZeroBias_pt) > 0)
events2=ak.sum(ak.num(Zmumu_pt) > 0)
#Prepare plot
data_2=[data1, data2]
bins=np.linspace(0.0,0.3,40)
colors=["#0072B2", "#FD0000"]
label1=f"Zero Bias, events={events1}"
label2=fr"Z $\to \mu \mu$, events={events2}"
labels=[label1, label2]
x_label=r"$\left( \frac{E_{iso}}{E_{t,muon}} \right)$"
y_label=r"Counts"
title = rf"$\left( \frac{{E_{{iso}}}}{{E_{{t,muon}}}}\right)$ ratio histogram, scaling={scaling}"
coolplot(data_2, bins, colors, labels, x_label, y_label, title, ax=axis[0, 1])

#------------------3rd subplot------------------------------------------------------------------------------------------------
#Set scaling factor
scaling=1.2
#Compute isolation
res1=muon_isolation_all_events(MuonTree_ZeroBias, ZeroBias_eta, ZeroBias_phi, dr_min, dr_max, [nmin,nmax], scaling=scaling)
res2=muon_isolation_all_events(MuonTree_Zmumu, Zmumu_eta, Zmumu_phi, dr_min, dr_max, [nmin,nmax], scaling=scaling)
#Compute ratio
data1=ak.flatten(res1)/ak.flatten(ZeroBias_pt)
data2=ak.flatten(res2)/ak.flatten(Zmumu_pt)
#Compute number of non empty events
events1=ak.sum(ak.num(ZeroBias_pt) > 0)
events2=ak.sum(ak.num(Zmumu_pt) > 0)
#Prepare plot
data_3=[data1, data2]
bins=np.linspace(0.0,0.3,40)
colors=["#0072B2", "#FD0000"]
label1=f"Zero Bias, events={events1}"
label2=fr"Z $\to \mu \mu$, events={events2}"
labels=[label1, label2]
x_label=r"$\left( \frac{E_{iso}}{E_{t,muon}} \right)$"
y_label=r"Counts"
title = rf"$\left( \frac{{E_{{iso}}}}{{E_{{t,muon}}}}\right)$ ratio histogram, scaling={scaling}"
coolplot(data_3, bins, colors, labels, x_label, y_label, title, ax=axis[1, 0])

#------------------4th subplot------------------------------------------------------------------------------------------------
#Set scaling factor
scaling=1.5
#Compute isolation
res1=muon_isolation_all_events(MuonTree_ZeroBias, ZeroBias_eta, ZeroBias_phi, dr_min, dr_max, [nmin,nmax], scaling=scaling)
res2=muon_isolation_all_events(MuonTree_Zmumu, Zmumu_eta, Zmumu_phi, dr_min, dr_max, [nmin,nmax], scaling=scaling)
#Compute ratio
data1=ak.flatten(res1)/ak.flatten(ZeroBias_pt)
data2=ak.flatten(res2)/ak.flatten(Zmumu_pt)
#Compute number of non empty events
events1=ak.sum(ak.num(ZeroBias_pt) > 0)
events2=ak.sum(ak.num(Zmumu_pt) > 0)
#Prepare plot
data_4=[data1, data2]
bins=np.linspace(0.0,0.3,40)
colors=["#0072B2", "#FD0000"]
label1=f"Zero Bias, events={events1}"
label2=fr"Z $\to \mu \mu$, events={events2}"
labels=[label1, label2]
x_label=r"$\left( \frac{E_{iso}}{E_{t,muon}} \right)$"
y_label=r"Counts"
title = rf"$\left( \frac{{E_{{iso}}}}{{E_{{t,muon}}}}\right)$ ratio histogram, scaling={scaling}"
coolplot(data_4, bins, colors, labels, x_label, y_label, title, ax=axis[1, 1])


#------------------5th subplot------------------------------------------------------------------------------------------------
#Set scaling factor
scaling=1.8
#Compute isolation
res1=muon_isolation_all_events(MuonTree_ZeroBias, ZeroBias_eta, ZeroBias_phi, dr_min, dr_max, [nmin,nmax], scaling=scaling)
res2=muon_isolation_all_events(MuonTree_Zmumu, Zmumu_eta, Zmumu_phi, dr_min, dr_max, [nmin,nmax], scaling=scaling)
#Compute ratio
data1=ak.flatten(res1)/ak.flatten(ZeroBias_pt)
data2=ak.flatten(res2)/ak.flatten(Zmumu_pt)
#Compute number of non empty events
events1=ak.sum(ak.num(ZeroBias_pt) > 0)
events2=ak.sum(ak.num(Zmumu_pt) > 0)
#Prepare plot
data_5=[data1, data2]
bins=np.linspace(0.0,0.3,40)
colors=["#0072B2", "#FD0000"]
label1=f"Zero Bias, events={events1}"
label2=fr"Z $\to \mu \mu$, events={events2}"
labels=[label1, label2]
x_label=r"$\left( \frac{E_{iso}}{E_{t,muon}} \right)$"
y_label=r"Counts"
title = rf"$\left( \frac{{E_{{iso}}}}{{E_{{t,muon}}}}\right)$ ratio histogram, scaling={scaling}"
coolplot(data_5, bins, colors, labels, x_label, y_label, title, ax=axis[2, 0])


#------------------6th subplot------------------------------------------------------------------------------------------------
#Set scaling factor
scaling=2
#Compute isolation
res1=muon_isolation_all_events(MuonTree_ZeroBias, ZeroBias_eta, ZeroBias_phi, dr_min, dr_max, [nmin,nmax], scaling=scaling)
res2=muon_isolation_all_events(MuonTree_Zmumu, Zmumu_eta, Zmumu_phi, dr_min, dr_max, [nmin,nmax], scaling=scaling)
#Compute ratio
data1=ak.flatten(res1)/ak.flatten(ZeroBias_pt)
data2=ak.flatten(res2)/ak.flatten(Zmumu_pt)
#Compute number of non empty events
events1=ak.sum(ak.num(ZeroBias_pt) > 0)
events2=ak.sum(ak.num(Zmumu_pt) > 0)
#Prepare plot
data_6=[data1, data2]
bins=np.linspace(0.0,0.3,40)
colors=["#0072B2", "#FD0000"]
label1=f"Zero Bias, events={events1}"
label2=fr"Z $\to \mu \mu$, events={events2}"
labels=[label1, label2]
x_label=r"$\left( \frac{E_{iso}}{E_{t,muon}} \right)$"
y_label=r"Counts"
title = rf"$\left( \frac{{E_{{iso}}}}{{E_{{t,muon}}}}\right)$ ratio histogram, scaling={scaling}"
coolplot(data_6, bins, colors, labels, x_label, y_label, title, ax=axis[2, 1])

plt.show()

# %%
fig, axis = plt.subplots(1,1 , figsize=(8, 6))
ROC_curve_compare_scaling(MuonTree_Zmumu,MuonTree_ZeroBias,Zmumu_pt,Zmumu_eta,Zmumu_phi,
                          ZeroBias_pt, ZeroBias_eta, ZeroBias_phi, event_range=[nmin,nmax], bin_range=[0,0.5], amount_of_curves=4)

# %%



