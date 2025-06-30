import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from my_functions import*

# %%
file= uproot.open("/home/ruben/Escritorio/BachelorThesisRuben/Data/Muon_trees.root") #opening the Root file with Uproot 

# %%
MuonTree_ZeroBias=file["MuonTree_ZeroBias;1"]
MuonTree_Zmumu=file["MuonTree_Zmumu;1"]

#Select range of events
nmin=0
nmax=10000

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
#Computo isolations
res1=muon_isolation_all_events(MuonTree_Zmumu, Zmumu_eta, Zmumu_phi, 0.09, 0.30, [nmin, nmax], int((nmax-nmin)/10))
res2=muon_isolation_all_events(MuonTree_ZeroBias, ZeroBias_eta, ZeroBias_phi, 0.09, 0.30, [nmin, nmax], int((nmax-nmin)/10))

# %%
#substract a global value of 5GeV to the isolation energies, and set negative results to 0
value=5*10**3

aux=ak.Array(res1)-value
res1_shift=ak.where(aux < 0, 0, aux)

aux=ak.Array(res2)-value
res2_shift=ak.where(aux < 0, 0, aux)

#Remove NaN
res1_shift=res1_shift[~np.isnan(res1_shift)]
res2_shift=res2_shift[~np.isnan(res2_shift)]

# %%
#Plot result
data1= ak.flatten(res1_shift)
data2= ak.flatten(res2_shift)
coolplot([data1, data2], np.linspace(0.0,3*10**4,20),
          labels=[fr"Z $\rightarrow \mu \mu$, events={len(data1)}", fr"ZeroBias, events={len(data2)}"], x_label="Isolation energy -5000 (MeV)",
          y_label="Counts", title=r"Isolation energy histogram ($E_{iso}$-5000 MeV)")

# %%
#Plot result
data1= ak.flatten(res1_shift)
data2= ak.flatten(res2_shift)
coolplot([data1[data1>0], data2[data2>0]], np.linspace(0.0,3*10**4,20),
          labels=[fr"Z $\rightarrow \mu \mu$, events={len(data1[data1>0])}", fr"ZeroBias, events={len(data2[data2>0])}"], x_label="Isolation energy -5000 (MeV)",
          y_label="Counts", title=r"Isolation energy histogram ($E_{iso}$-5000 MeV) - Only > 0 values")

# %%
#Compute ratio
ratio1=ak.flatten(res1_shift)/ak.flatten(Zmumu_pt[~np.isnan(res1_shift)])
ratio2=ak.flatten(res2_shift)/ak.flatten(ZeroBias_pt[~np.isnan(res2_shift)])

# %%
#Plot result
data1= ratio1
data2= ratio2
coolplot([ratio1, ratio2], np.linspace(0.0,0.5,20),labels=[fr"Z $\rightarrow \mu \mu$, events={len(ratio1)}",
 fr"ZeroBias, events={len(ratio2)}"], x_label=r"$\left(\frac{E_{iso} -5000 (MeV)}{E_{t}}\right)$ ratio", y_label="Counts", 
 title=r"$\left(\frac{E_{iso} -5000 (MeV)}{E_{t}}\right)$ ratio histogram")



# %%
#Plot result
data1= ratio1
data2= ratio2
coolplot([ratio1[ratio1>0], ratio2[ratio2>0]], np.linspace(0.0,0.5,20),labels=[fr"Z $\rightarrow \mu \mu$, muons={len(ratio1[ratio1>0])}",
 fr"ZeroBias, muons={len(ratio2[ratio2>0])}"], x_label=r"$\left(\frac{E_{iso} -5000 (MeV)}{E_{t}}\right)$ ratio", y_label="Counts", 
 title=r"$\left(\frac{E_{iso} -5000 (MeV)}{E_{t}}\right)$ ratio histogram - Only > 0 values")
