import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from my_functions import *

file= uproot.open("/home/ruben/Escritorio/BachelorThesisRuben/Data/Muon_trees.root") #opening the Root file with Uproot 
MuonTree_ZeroBias=file["MuonTree_ZeroBias;1"]
MuonTree_Zmumu=file["MuonTree_Zmumu;1"]

# %%
#Prepare the data for the plots

#Choose the range of events to plot
nmin=0
nmax=10000
dr_min=0.09
dr_max=0.3

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
res=muon_isolation_all_events(MuonTree_ZeroBias,ZeroBias_eta,ZeroBias_phi,dr_min,dr_max,[nmin,nmax],1500)
data1=ak.flatten(res)

#Compute the isolation and prepare it for plotting Z mu mu
res=muon_isolation_all_events(MuonTree_Zmumu,Zmumu_eta,Zmumu_phi,dr_min,dr_max,[nmin,nmax],1500)
data2=ak.flatten(res)

# %%
#Prepare the data
ratio1=data1/ak.flatten(ZeroBias_pt)
ratio2=data2/ak.flatten(Zmumu_pt)

e1=ak.flatten(ZeroBias_pt)
e2=ak.flatten(Zmumu_pt)
#Prepare the limits
xlim=100000
ylim=0.35

#Remove NaN and empty
mask1 = ~np.isnan(data1) & ~np.isnan(ratio1) & ~np.isinf(data1) & ~np.isinf(ratio1)
e1= e1[mask1]
ratio1= ratio1[mask1]
mask2 = ~np.isnan(data2) & ~np.isnan(ratio2) & ~np.isinf(data2) & ~np.isinf(ratio2)
e2= e2[mask2]
ratio2= ratio2[mask2]

#Perform a linear regression
slope1, intercept1 = np.polyfit(e1,ratio1,1)
slope2, intercept2 = np.polyfit(e2,ratio2,1)
x_regression=np.linspace(0,xlim,100)
y_regression1=slope1*x_regression+intercept1
y_regression2=slope2*x_regression+intercept2


# Convert Awkward arrays to NumPy
e1_np = ak.to_numpy(e1)
ratio1_np = ak.to_numpy(ratio1)
e2_np = ak.to_numpy(e2)
ratio2_np = ak.to_numpy(ratio2)

#Plot heatmaps
# Set up the figure and axes
fig, axs = plt.subplots(1, 2, figsize=(14, 6), sharey=True)

bins = 80

# --- Zero Bias heatmap ---
h1 = axs[0].hist2d(e1_np, ratio1_np, bins=bins, range=[[10000, xlim], [0, ylim]],
                   cmap="Blues", norm=plt.cm.colors.LogNorm())
axs[0].plot(x_regression, y_regression1, color="black",
            label=fr"Zero Bias: y = {slope1:.2e}x + {intercept1:.2f}")
axs[0].set_xlabel("Transverse Energy (MeV)")
axs[0].set_ylabel(r"$\frac{E_{T,iso}}{E_{T}}$")
axs[0].set_title("Zero Bias")
axs[0].legend()
axs[0].grid(alpha=0.5, linestyle="--")
axs[0].set_xlim(10000, xlim)
axs[0].set_ylim(0, ylim)
fig.colorbar(h1[3], ax=axs[0], label='Counts')

# --- Z → μμ heatmap ---
h2 = axs[1].hist2d(e2_np, ratio2_np, bins=bins, range=[[10000, xlim], [0, ylim]],
                   cmap="Reds", norm=plt.cm.colors.LogNorm())
axs[1].plot(x_regression, y_regression2, color="black",
            label=fr"Z $\to \mu\mu$: y = {slope2:.2e}x + {intercept2:.2f}")
axs[1].set_xlabel("Transverse Energy (MeV)")
axs[1].set_title(r"Z $\to \mu\mu$")
axs[1].legend()
axs[1].grid(alpha=0.5, linestyle="--")
axs[1].set_xlim(10000, xlim)
axs[1].set_ylim(0, ylim)
fig.colorbar(h2[3], ax=axs[1], label='Counts')

# Adjust layout
fig.suptitle(r"Heatmaps: $\frac{E_{T,iso}}{E_T}$ vs $E_T$", fontsize=14)
plt.tight_layout(rect=[0, 0.03, 1, 1])
plt.show()