import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
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
Zmumu_e=quality_selector_with_empty(MuonTree_Zmumu["muon_quality"].array(),MuonTree_Zmumu["muon_e"].array(),0)[nmin1:nmax1]
Zmumu_charge=quality_selector_with_empty(MuonTree_Zmumu["muon_quality"].array(),MuonTree_Zmumu["muon_charge"].array(),0)[nmin1:nmax1]
#Select the Z peak pairs
Zmumu_pt, Zmumu_eta, Zmumu_phi= get_all_Z_peak_pairs(Zmumu_pt,Zmumu_eta,Zmumu_phi)

#Open ZeroBias file
file= uproot.open("/home/ruben/Escritorio/BachelorThesisRuben/Data/zbV3_skim.root") 
MuonTree_ZeroBias=file["MuonTree;1"]

#Apply energy cut to offline
ZeroBias_pt=energy_cut_with_empty(MuonTree_ZeroBias["muon_pt"].array(), MuonTree_ZeroBias["muon_pt"].array())
ZeroBias_eta=energy_cut_with_empty(MuonTree_ZeroBias["muon_pt"].array(), MuonTree_ZeroBias["muon_eta"].array())
ZeroBias_phi=energy_cut_with_empty(MuonTree_ZeroBias["muon_pt"].array(), MuonTree_ZeroBias["muon_phi"].array())
ZeroBias_e=energy_cut_with_empty(MuonTree_ZeroBias["muon_pt"].array(), MuonTree_ZeroBias["muon_e"].array())
ZeroBias_charge=energy_cut_with_empty(MuonTree_ZeroBias["muon_pt"].array(), MuonTree_ZeroBias["muon_charge"].array())

#Get online data
Zbl1_pt=MuonTree_ZeroBias["LVL1Muon_et"].array() * 1000
Zbl1_eta=MuonTree_ZeroBias["LVL1Muon_eta"].array()
Zbl1_phi=MuonTree_ZeroBias["LVL1Muon_phi"].array()

#Create mask matching offline and LVL1
mask=offline_LVL1_matcher(ZeroBias_eta, ZeroBias_phi, Zbl1_eta, Zbl1_phi)

#Apply mask
ZeroBias_pt=ZeroBias_pt[mask]
ZeroBias_eta=ZeroBias_eta[mask]
ZeroBias_phi=ZeroBias_phi[mask]
ZeroBias_e=ZeroBias_e[mask]
ZeroBias_charge=ZeroBias_charge[mask]

data1=Zmumu_pt
data2=ZeroBias_pt
l1=len(ak.flatten(data1))
l2=len(ak.flatten(data2))
data=[data1,data2]
bins=np.linspace(0,1*10**5,50)
colors=['r','b']
x_label=r"Muon transverse momentum $p_T$ (MeV)"
y_label="Counts"
title=r"$p_T$ histogram"
label1=rf'Z $\longrightarrow \mu \mu$, muons={l1}'
label2=f'Zero Bias, muons={l2}'
labels=[label1,label2]

coolplot(data,bins,colors,labels,x_label,y_label,title, collect_overflow=False, plot_show=False)
plt.savefig('pt_hist.pdf', format='pdf')
plt.show()
# %%
data1=Zmumu_eta
data2=ZeroBias_eta
l1=len(ak.flatten(data1))
l2=len(ak.flatten(data2))

data=[data1,data2]
bins=np.linspace(-4,4,50)
colors=['r','b']
x_label=r"Muon pseudorapitidy $\eta$"
y_label="Counts"
title=r"Muon counts vs $\eta$ ATLAS Detector"
label1=rf'Z $\longrightarrow \mu \mu$ data, muons={l1}'
label2=f'Zero Bias data, muons={l2}'
labels=[label1,label2]

coolplot(data,bins,colors,labels,x_label,y_label,title, collect_overflow=False, plot_show=False)
plt.savefig('eta_hist.pdf', format='pdf')
plt.show()
# %%
data1=Zmumu_phi
data2=ZeroBias_phi
l1=len(ak.flatten(data1))
l2=len(ak.flatten(data2))

data=[data1,data2]
bins=np.linspace(-4,4,50)
colors=['r','b']
x_label=r"Muon phi $\phi$"
y_label="Counts"
title=r"Muon counts vs $\phi$ ATLAS Detector"
label1=rf'Z $\longrightarrow \mu \mu$ data, muons={l1}'
label2=f'Zero Bias data, muons={l2}'
labels=[label1,label2]

coolplot(data,bins,colors,labels,x_label,y_label,title, collect_overflow=False, plot_show=False)
plt.savefig('phi_hist.pdf', format='pdf')
plt.show()
# %%
data1=Zmumu_e
data2=ZeroBias_e
l1=len(ak.flatten(data1))
l2=len(ak.flatten(data2))

data=[data1,data2]
bins=np.linspace(0,2*10**5,50)
colors=['r','b']
x_label=r"Muon Energy $E$ (MeV)"
y_label="Counts"
title=r"Muon counts vs $E$ ATLAS Detector"
label1=rf'Z $\longrightarrow \mu \mu$ data, muons={l1}'
label2=f'Zero Bias data, muons={l2}'
labels=[label1,label2]

coolplot(data,bins,colors,labels,x_label,y_label,title, collect_overflow=False, plot_show=False)
plt.savefig('e_hist.pdf', format='pdf')
plt.show()
# %%
data1=Zmumu_charge
data2=ZeroBias_charge
l1=len(ak.flatten(data1))
l2=len(ak.flatten(data2))

data=[data1,data2]
bins=np.linspace(-2,2,50)
colors=['r','b']
x_label=r"Muon charge (e-)"
y_label="Counts"
title=r"Muon counts vs charge ATLAS Detector"
label1=rf'Z $\longrightarrow \mu \mu$ data, muons={l1}'
label2=f'Zero Bias data, muons={l2}'
labels=[label1,label2]

coolplot(data,bins,colors,labels,x_label,y_label,title, collect_overflow=False, plot_show=False)
plt.savefig('charge_hist.pdf', format='pdf')
plt.show()