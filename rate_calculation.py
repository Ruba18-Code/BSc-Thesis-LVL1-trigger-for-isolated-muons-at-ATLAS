# %%
from my_functions import *

# %%
file= uproot.open("/home/ruben/Escritorio/BachelorThesisRuben/Data/Muon_trees.root") #opening the Root file with Uproot 

MuonTree_Zmumu=file["MuonTree_Zmumu;1"]

Zmumu_data=MuonTree_Zmumu["muon_e"].array()

MuonTree_Zmumu=file["MuonTree_Zmumu;1"] 
MuonTree_ZeroBias=file["MuonTree_ZeroBias;1"]


et=MuonTree_ZeroBias["LVL1Muon_et"].array()

# %%
scaling_factor=(2340/3564)*40000 #Scaling factor in kHz!!!
x=np.linspace(0,18,90)

res=[]
for val in x:
    a, _ = rate_calculator(et,val,scaling_factor)
    res.append(a)

plt.plot(x,res,color='b')
plt.tight_layout()
plt.legend(["LVL1 Zero Bias"])
plt.grid(True, linestyle='--', alpha=0.5)
plt.xlabel(r"Muon $e_T$ (MeV)")
plt.ylabel("Rate (kHz)")
plt.title(r"LVL1 Muon rate vs $e_T$")
plt.show()

# %%
scaling_factor=(2340/3564)*40000  #in kHz
val=14000
x, y = rate_calculator(et,14,scaling_factor)
expected_rate=10 #kHz

frac=x/expected_rate

print("The estimated rate for a threshold of", val, "MeV is", x,"kHz. This is", frac, "times bigger than the expected value,", expected_rate, "kHz")

#For the Zero bias data we do not select quality (because it's noise anyways)

data=[et]
bins=np.linspace(0,18,18)
colors=['b']
x_label=r"Muon transverse energy ($e_T$)"
y_label="Counts"
title=r"LVL1 Muon rate vs $e_T$ ATLAS Detector"
label2='LVL1 Zero Bias'
labels=[label2]

coolplot(data,bins,colors,labels,x_label,y_label,title)