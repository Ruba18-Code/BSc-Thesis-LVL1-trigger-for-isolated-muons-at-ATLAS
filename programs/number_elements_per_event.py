import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

#Import functions------------------------------------------------------------------------------------------------------
from my_functions import*
#Open ROOT file with uproot--------------------------------------------------------------------------------------------
file= uproot.open("/home/ruben/Escritorio/BachelorThesisRuben/Data/Muon_trees.root") 
#Get trees-------------------------------------------------------------------------------------------------------------
MuonTree_ZeroBias=file["MuonTree_ZeroBias;1"]
MuonTree_Zmumu=file["MuonTree_Zmumu;1"]
#Get data from trees---------------------------------------------------------------------------------------------------
#Z -> mu mu
Zmumu_pt  = MuonTree_Zmumu["muon_pt"].array()
Zmumu_eta = MuonTree_Zmumu["muon_eta"].array()
Zmumu_phi = MuonTree_Zmumu["muon_phi"].array()
#Zero Bias
ZeroBias_pt  = MuonTree_ZeroBias["muon_pt"].array()
ZeroBias_eta = MuonTree_ZeroBias["muon_eta"].array()
ZeroBias_phi = MuonTree_ZeroBias["muon_phi"].array()
#Pre-select data-------------------------------------------------------------------------------------------------------
#Chose event range
#Z -> mu mu range
nmin1=0
nmax1=50000
#Select quality 0 Z->mumu
Zmumu_pt=quality_selector_with_empty(MuonTree_Zmumu["muon_quality"].array(),MuonTree_Zmumu["muon_pt"].array(),0)[nmin1:nmax1]
Zmumu_eta=quality_selector_with_empty(MuonTree_Zmumu["muon_quality"].array(),MuonTree_Zmumu["muon_eta"].array(),0)[nmin1:nmax1]
Zmumu_phi=quality_selector_with_empty(MuonTree_Zmumu["muon_quality"].array(),MuonTree_Zmumu["muon_phi"].array(),0)[nmin1:nmax1]

Zmumu_num=ak.num(Zmumu_pt)
#Pre-select data-------------------------------------------------------------------------------------------------------
Zmumu_pt=quality_selector_with_empty(MuonTree_Zmumu["muon_quality"].array(),MuonTree_Zmumu["muon_pt"].array(),0)[nmin1:nmax1]
Zmumu_eta=quality_selector_with_empty(MuonTree_Zmumu["muon_quality"].array(),MuonTree_Zmumu["muon_eta"].array(),0)[nmin1:nmax1]
Zmumu_phi=quality_selector_with_empty(MuonTree_Zmumu["muon_quality"].array(),MuonTree_Zmumu["muon_phi"].array(),0)[nmin1:nmax1]
#Select the Z peak pairs
Zmumu_pt, Zmumu_eta, Zmumu_phi= get_all_Z_peak_pairs(Zmumu_pt,Zmumu_eta,Zmumu_phi)
# Determine the length for the Zmumu histogram label
Zmumu_num_after=ak.num(Zmumu_pt)
l2 = len(ak.flatten(Zmumu_pt))
Zmumu_num_after=Zmumu_num_after[Zmumu_num_after > 0]

l1 = len(Zmumu_num) 
l2 = len(Zmumu_num_after)  

# Set up the figure and subplots
fig, axs = plt.subplots(1, 2, figsize=(12, 6))  

# Define bins
bins = np.arange(0, 10, 1)

# Plot Zmumu histogram on the first subplot
axs[0].hist(Zmumu_num, bins=bins, edgecolor='black', color='r')
axs[0].set_ylabel("Counts")
axs[0].set_xlabel("Number of particles per event")
axs[0].legend()
axs[0].set_title(rf"Z$\to\mu\mu$, n={l1}")
axs[0].grid(alpha=0.5, linestyle='--')

# Plot Zmumu histogram after pre-selection on the second subplot
axs[1].hist(Zmumu_num_after, bins=bins, edgecolor='black', color='r')
axs[1].set_ylabel("Counts")
axs[1].set_xlabel("Number of particles per event")
axs[1].legend()
axs[1].set_title(rf"Z$\to\mu\mu$ after pre-selection, n={l2}")
axs[1].grid(alpha=0.5, linestyle='--')

# Adjust layout to avoid overlap
plt.tight_layout()
plt.savefig('plot.pdf', format='pdf')
# Display the plot
plt.show()
