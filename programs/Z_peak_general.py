import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from my_functions import*


file= uproot.open("/home/ruben/Escritorio/BachelorThesisRuben/Data/Muon_trees.root") #opening the Root file with Uproot 
MuonTree_Zmumu=file["MuonTree_Zmumu;1"]
#Z -> mu mu range
nmin1=0
nmax1=200000


Zmumu_pt=MuonTree_Zmumu["muon_pt"].array()[nmin1:nmax1]
Zmumu_eta=MuonTree_Zmumu["muon_eta"].array()[nmin1:nmax1]
Zmumu_phi=MuonTree_Zmumu["muon_phi"].array()[nmin1:nmax1]

Z_mass1=invariant_mass_all_muons(Zmumu_pt,Zmumu_eta,Zmumu_phi) #computes the invariant mass vector
Z_mass1=Z_mass1[~np.isnan(Z_mass1)]
n1=len(Z_mass1)
#Select quality 0 Z->mumu
Zmumu_pt=quality_selector(MuonTree_Zmumu["muon_quality"].array(),MuonTree_Zmumu["muon_pt"].array(),0)[nmin1:nmax1]
Zmumu_eta=quality_selector(MuonTree_Zmumu["muon_quality"].array(),MuonTree_Zmumu["muon_eta"].array(),0)[nmin1:nmax1]
Zmumu_phi=quality_selector(MuonTree_Zmumu["muon_quality"].array(),MuonTree_Zmumu["muon_phi"].array(),0)[nmin1:nmax1]

Z_mass2=invariant_mass_all_muons(Zmumu_pt,Zmumu_eta,Zmumu_phi) #computes the invariant mass vector
Z_mass2=Z_mass2[~np.isnan(Z_mass2)]
n2=len(Z_mass2)
#COmpute standard deviations
mean1=ak.mean(Z_mass1)
std1 = np.std(Z_mass1, ddof=1)

mean2=ak.mean(Z_mass2)
std2 = np.std(Z_mass2, ddof=1)

data=[Z_mass1, Z_mass2]
bins=np.linspace(6*10**4,1.2*10**5, 40)
colors=['r','b']
x_label=r"Energy (MeV)"
y_label="Counts"
title=rf"Invariant mass histogram"
label1=f'Any quality - mean={int(mean1)},\n std={int(std1)}, n={n1}'
label2=f"Quality=0 - mean={int(mean2)},\n std={int(std2)}, n={n2}"
labels=[label1, label2]
coolplot(data,bins,colors,labels,x_label,y_label,title)
plt.savefig('Zpeak.pdf', format='pdf')