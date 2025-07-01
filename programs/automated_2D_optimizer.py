
# %%
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from my_functions import*

# %%
file= uproot.open("/home/ruben/Escritorio/BachelorThesisRuben/Data/Muon_trees.root") #opening the Root file with Uproot 

MuonTree_ZeroBias=file["MuonTree_ZeroBias;1"]
MuonTree_Zmumu=file["MuonTree_Zmumu;1"]

#Select range of events
nmin1=0
nmax1=4000

#Select quality 0 Z->mumu
Zmumu_pt=quality_selector_with_empty(MuonTree_Zmumu["muon_quality"].array(),MuonTree_Zmumu["muon_pt"].array(),0)[nmin1:nmax1]
Zmumu_eta=quality_selector_with_empty(MuonTree_Zmumu["muon_quality"].array(),MuonTree_Zmumu["muon_eta"].array(),0)[nmin1:nmax1]
Zmumu_phi=quality_selector_with_empty(MuonTree_Zmumu["muon_quality"].array(),MuonTree_Zmumu["muon_phi"].array(),0)[nmin1:nmax1]
#And select the Z peak pairs
Zmumu_pt1, Zmumu_eta1, Zmumu_phi1= get_all_Z_peak_pairs(Zmumu_pt,Zmumu_eta,Zmumu_phi)

nmin2=0
nmax2=100000
ZeroBias_eta1=energy_cut_with_empty(MuonTree_ZeroBias["muon_pt"].array(), MuonTree_ZeroBias["muon_eta"].array())[nmin2:nmax2]
ZeroBias_phi1=energy_cut_with_empty(MuonTree_ZeroBias["muon_pt"].array(), MuonTree_ZeroBias["muon_phi"].array())[nmin2:nmax2]
ZeroBias_pt1=energy_cut_with_empty(MuonTree_ZeroBias["muon_pt"].array(), MuonTree_ZeroBias["muon_pt"].array())[nmin2:nmax2]

##################################################################################################################################33
points=1
iterations=1
next_dr_mins=np.linspace(0,0.2,points)
next_dr_maxs=np.linspace(0.3,1,points)

min_range=[min(next_dr_mins), max(next_dr_mins)]
max_range=[min(next_dr_maxs), max(next_dr_maxs)]

FPR_effs, dr_mins, dr_maxs= ROC_FPR_2D_plot(MuonTree_Zmumu, MuonTree_ZeroBias, Zmumu_pt1, Zmumu_eta1, Zmumu_phi1, ZeroBias_pt1, ZeroBias_eta1,
                 ZeroBias_phi1, [nmin1,nmax1], [nmin2, nmax2], min_range, max_range, points)

#Get flat position the first 5 indices of the sorted array
flat_indices = np.argsort(FPR_effs, axis=None)[:5] 
#Arrange them into 2D coordinates
positions = np.unravel_index(flat_indices, FPR_effs.shape) 
#Create list of respective pairs 
best_coords = list(zip(positions[0], positions[1]))
#Print and prepare next iteration
next_dr_mins=[]
next_dr_maxs=[]
print("Top 5 lowest FPR(90%):")
for i, (row, col) in enumerate(best_coords):
    print(fr"{i+1}. ΔR = [{dr_mins[row]}, {dr_maxs[col]}] → FPR = {FPR_effs[row, col]}")
    next_dr_mins.append(dr_mins[row])
    next_dr_maxs.append(dr_maxs[col])

for i in range(iterations):
    min_range=[min(next_dr_mins), max(next_dr_mins)]
    max_range=[min(next_dr_maxs), max(next_dr_maxs)]

    FPR_effs, dr_mins, dr_maxs= ROC_FPR_2D_plot(MuonTree_Zmumu, MuonTree_ZeroBias, Zmumu_pt1, Zmumu_eta1, Zmumu_phi1, ZeroBias_pt1, ZeroBias_eta1,
                    ZeroBias_phi1, [nmin1,nmax1], [nmin2, nmax2], min_range, max_range, points)

    #Get flat position the first 5 indices of the sorted array
    flat_indices = np.argsort(FPR_effs, axis=None)[:5] 
    #Arrange them into 2D coordinates
    positions = np.unravel_index(flat_indices, FPR_effs.shape) 
    #Create list of respective pairs 
    best_coords = list(zip(positions[0], positions[1]))
    #Print and prepare next iteration
    next_dr_mins=[]
    next_dr_maxs=[]
    print("Top 5 lowest FPR(90%):")
    for i, (row, col) in enumerate(best_coords):
        print(fr"{i+1}. ΔR = [{dr_mins[row]}, {dr_maxs[col]}] → FPR = {FPR_effs[row, col]}")
        next_dr_mins.append(dr_mins[row])
        next_dr_maxs.append(dr_maxs[col])