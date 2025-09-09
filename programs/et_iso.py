import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
#Import functions------------------------------------------------------------------------------------------------------
from my_functions import*
#Open Zmumu file
#Import functions------------------------------------------------------------------------------------------------------
from my_functions import*
#Open Zmumu file
file= uproot.open("/home/ruben/Escritorio/BachelorThesisRuben/Data/Muon_trees.root") 
MuonTree_Zmumu=file["MuonTree_Zmumu;1"]

#Set event range
nmin1=0
nmax1=3000

#Choose quality 0
Zmumu_pt=quality_selector_with_empty(MuonTree_Zmumu["muon_quality"].array(),MuonTree_Zmumu["muon_pt"].array(),0)[nmin1:nmax1]
Zmumu_eta=quality_selector_with_empty(MuonTree_Zmumu["muon_quality"].array(),MuonTree_Zmumu["muon_eta"].array(),0)[nmin1:nmax1]
Zmumu_phi=quality_selector_with_empty(MuonTree_Zmumu["muon_quality"].array(),MuonTree_Zmumu["muon_phi"].array(),0)[nmin1:nmax1]
#Select the Z peak pairs
Zmumu_pt, Zmumu_eta, Zmumu_phi= get_all_Z_peak_pairs(Zmumu_pt,Zmumu_eta,Zmumu_phi)

#Open ZeroBias file
file= uproot.open("/home/ruben/Escritorio/BachelorThesisRuben/Data/zbV3_skim.root") 
MuonTree_ZeroBias=file["MuonTree;1"]

#Apply energy cut to offline
ZeroBias_pt=energy_cut_with_empty(MuonTree_ZeroBias["muon_pt"].array(), MuonTree_ZeroBias["muon_pt"].array())
ZeroBias_eta=energy_cut_with_empty(MuonTree_ZeroBias["muon_pt"].array(), MuonTree_ZeroBias["muon_eta"].array())
ZeroBias_phi=energy_cut_with_empty(MuonTree_ZeroBias["muon_pt"].array(), MuonTree_ZeroBias["muon_phi"].array())

#Get online data
Zbl1_pt=MuonTree_ZeroBias["LVL1Muon_et"].array() * 1000
Zbl1_eta=MuonTree_ZeroBias["LVL1Muon_eta"].array()
Zbl1_phi=MuonTree_ZeroBias["LVL1Muon_phi"].array()

#Create mask matching offline and LVL1
mask=offline_LVL1_matcher(ZeroBias_eta, ZeroBias_phi, Zbl1_eta, Zbl1_phi, dr_threshold=0.4)

#Apply mask
ZeroBias_pt=ZeroBias_pt[mask]
ZeroBias_eta=ZeroBias_eta[mask]
ZeroBias_phi=ZeroBias_phi[mask]

dr_min=0.1
dr_max=0.3
#Compute the isolation and prepare it for plotting ZeroBias
res=muon_isolation_all_events(MuonTree_ZeroBias,ZeroBias_eta,ZeroBias_phi,dr_min,dr_max,[0, len(ZeroBias_pt)])
data1=ak.flatten(res)

#Compute the isolation and prepare it for plotting Z mu mu
res=muon_isolation_all_events(MuonTree_Zmumu,Zmumu_eta,Zmumu_phi,dr_min,dr_max,[nmin1,nmax1])
data2=ak.flatten(res)

# %%
#Prepare the data
ratio1=data1
ratio2=data2

e1=ak.flatten(ZeroBias_pt)
e2=ak.flatten(Zmumu_pt)
#Prepare the limits
xlim=100000
ylim=15000

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

axs[0].set_xlabel("Transverse Energy [MeV]")
axs[0].set_ylabel(r"$E_{T,iso}$")
axs[0].set_title(f"Zero Bias, muons={len(e1_np)}")
axs[0].legend()
axs[0].grid(alpha=0.5, linestyle="--")
axs[0].set_xlim(10000, xlim)
axs[0].set_ylim(0, ylim)
fig.colorbar(h1[3], ax=axs[0], label='Counts')

# --- Z → μμ heatmap ---
h2 = axs[1].hist2d(e2_np, ratio2_np, bins=bins, range=[[10000, xlim], [0, ylim]],
                   cmap="Reds", norm=plt.cm.colors.LogNorm())

axs[1].set_xlabel("Transverse Energy [MeV]")
axs[1].set_ylabel(r"$E_{T,iso}$")
axs[1].set_title(rf"Z $\to \mu\mu$, muons={len(e2_np)}")
axs[1].legend()
axs[1].grid(alpha=0.5, linestyle="--")
axs[1].set_xlim(10000, xlim)
axs[1].set_ylim(0, ylim)
fig.colorbar(h2[3], ax=axs[1], label='Counts')

# Adjust layout
fig.suptitle(rf"Heatmaps: $E_{{T,iso}}$ vs $E_T$, $\Delta R$=[{dr_min},{dr_max}]", fontsize=14)
plt.tight_layout(rect=[0, 0.03, 1, 1])
plt.savefig('heatmaps_et_iso', format='pdf')
plt.show()
