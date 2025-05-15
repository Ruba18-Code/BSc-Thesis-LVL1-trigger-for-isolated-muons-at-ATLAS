# %%
from my_functions import *

# %%
file= uproot.open("/home/ruben/Escritorio/BachelorThesisRuben/Data/Muon_trees.root") #opening the Root file with Uproot 

MuonTree_Zmumu=file["MuonTree_Zmumu;1"]

Zmumu_data=MuonTree_Zmumu["muon_e"].array()

MuonTree_Zmumu=file["MuonTree_Zmumu;1"] 
MuonTree_ZeroBias=file["MuonTree_ZeroBias;1"]

eta=MuonTree_Zmumu["muon_eta"].array()
phi=MuonTree_Zmumu["muon_phi"].array()
pt=MuonTree_ZeroBias["muon_pt"].array()

# %%
scaling_factor=(2340/3564)*40000 #Scaling factor in kHz!!!
x=np.linspace(0,10**5,100)

res=[]
for val in x:
    a, _ = rate_calculator(pt,val,scaling_factor)
    res.append(a)

plt.plot(x,res,color='b')
plt.tight_layout()
plt.legend(["Zero Bias"])
plt.grid(True, linestyle='--', alpha=0.5)
plt.xlabel(r"Muon $p_T$ (MeV)")
plt.ylabel("Rate (kHz)")
plt.title(r"Muon rate vs $p_T$")
plt.show()

# %%
scaling_factor=(2340/3564)*40000  #in kHz
val=20000
x, _ = rate_calculator(pt,20000,scaling_factor)
expected_rate=10 #kHz

frac=x/expected_rate

print("The estimated rate for a threshold of", val, "MeV is", x,"kHz. This is", frac, "times bigger than the expected value,", expected_rate, "kHz")



