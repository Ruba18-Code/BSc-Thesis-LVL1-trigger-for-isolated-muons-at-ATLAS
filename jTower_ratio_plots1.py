# %%
from my_functions import*

# %%
file= uproot.open("/home/ruben/Escritorio/BachelorThesisRuben/Data/Muon_trees.root") #opening the Root file with Uproot 
MuonTree_ZeroBias=file["MuonTree_ZeroBias;1"]
MuonTree_Zmumu=file["MuonTree_Zmumu;1"]

# %%
#Set delta R values and number of events to plot

dr_min=0.0
dr_max=0.4

nmin1=0
nmax1=1000

nmin2=0
nmax2=1000

ZeroBias_eta=energy_cut(MuonTree_ZeroBias["muon_pt"].array(), MuonTree_ZeroBias["muon_eta"].array())[nmin2:nmax2]
ZeroBias_phi=energy_cut(MuonTree_ZeroBias["muon_pt"].array(), MuonTree_ZeroBias["muon_phi"].array())[nmin2:nmax2]
ZeroBias_pt=energy_cut(MuonTree_ZeroBias["muon_pt"].array(), MuonTree_ZeroBias["muon_pt"].array())[nmin2:nmax2]

# %%
"""
Z->mumu any quality
"""



#Compute the isolation for Z->mu mu events
pt_events=MuonTree_Zmumu["muon_pt"].array()[nmin1:nmax1]
eta_events=MuonTree_Zmumu["muon_eta"].array()[nmin1:nmax1]
phi_events=MuonTree_Zmumu["muon_phi"].array()[nmin1:nmax1]

data1=ak.flatten(pt_events)
isolation =muon_isolation_all_events(MuonTree_Zmumu, eta_events, phi_events, dr_min, dr_max, [nmin1,nmax1], 300)
data2=ak.flatten(isolation)

ratio1=data2/data1

#Get the number of events with muons
not_empty_count1 = ak.sum(ak.num(pt_events) > 0)

#Compute the isolation for Zero Bias events



data1=ak.flatten(ZeroBias_pt)
isolation =muon_isolation_all_events(MuonTree_ZeroBias, ZeroBias_eta, ZeroBias_phi, dr_min, dr_max, [nmin2,nmax2], 1000)
data2=ak.flatten(isolation)

ratio2=data2/data1

#Get the number of events with muons
not_empty_count2 = ak.sum(ak.num(pt_events) > 0)


"""
Z -> mu mu vs ZeroBias offline QUALITY 0
"""

#Compute the isolation for Z->mu mu events
pt_events=quality_selector(MuonTree_Zmumu["muon_quality"].array(), MuonTree_Zmumu["muon_pt"].array(), 0)[nmin1:nmax1]
eta_events=quality_selector(MuonTree_Zmumu["muon_quality"].array(), MuonTree_Zmumu["muon_eta"].array(), 0)[nmin1:nmax1]
phi_events=quality_selector(MuonTree_Zmumu["muon_quality"].array(), MuonTree_Zmumu["muon_phi"].array(), 0)[nmin1:nmax1]

data1=ak.flatten(pt_events)
isolation =muon_isolation_all_events(MuonTree_Zmumu, eta_events, phi_events, dr_min, dr_max, [nmin1,nmax1], 300)
data2=ak.flatten(isolation)

ratio3=data2/data1

#Get the number of events with muons
not_empty_count3 = ak.sum(ak.num(pt_events) > 0)

#Since we are using the same data, we can just copy the ratio and the not_empty_count
ratio4=ratio2
not_empty_count4 =not_empty_count2



# %%
"""
I define this functions to plot everything more comfortably and as subplots
"""

def f(i,ax):
    #This line is used to set the axis to the current axis
    plt.sca(ax)
    
    if i==0:
        #Set the colors and labels for the plot
        colors=["#0072B2", "#FD0000"]
        labels=[fr"Z $\longrightarrow \mu \mu$ (ANY Q), events={not_empty_count1}",fr"Zero Bias, events={not_empty_count2}"]

        #Plot the data
        coolplot([ratio1,ratio2],
                        np.linspace(0,6,20),
                        colors,labels,
                        "Isolation / Transverse energy","Counts",
                "ANY quality case",
                plot_show=False)

    if i==1:
        #Set the colors and labels for the plot
        colors=["#0072B2", "#FD0000"]
        labels=[fr"Z $\longrightarrow \mu \mu$ (Q=0), events={not_empty_count3}",fr"Zero Bias, events={not_empty_count4}"]

        #Plot the data
        coolplot([ratio3,ratio4],
                        np.linspace(0,6,20),
                        colors,labels,
                        "Isolation / Transverse energy","Counts",
                "Quality 0 case",
                plot_show=False)

def f_subplots():
    # Plot as subplots
    fig, ax = plt.subplots(1, 2, figsize=(18, 7))
    
    ax=ax.flatten()
    #for each lower and upper limit, plot the data
    for i in [0,1]:
        f(i,ax[i])

    fig.suptitle(fr"$\left(\frac{{E_{{t,\;iso}}}}{{E_t}}\right)$ ratio histogram, $\Delta R$ = [{dr_min}, {dr_max}]",fontsize=16)  # Set the global title
    plt.tight_layout(rect=[0, 0, 1, 0.98])  # Adjust layout to leave space for suptitle
    plt.show()

f_subplots()

# %%
"""
Z->mumu any quality
"""
nmin1=0
nmax1=1000
#Compute the isolation for Z->mu mu events
pt_events=MuonTree_Zmumu["LVL1Muon_et"].array()[nmin1:nmax1]
eta_events=MuonTree_Zmumu["LVL1Muon_eta"].array()[nmin1:nmax1]
phi_events=MuonTree_Zmumu["LVL1Muon_phi"].array()[nmin1:nmax1]

data1=ak.flatten(pt_events)*1000 #multiply by 1000 to convert to MeV
isolation =muon_isolation_all_events(MuonTree_Zmumu, eta_events, phi_events, dr_min, dr_max, [nmin1,nmax1], 300)
data2=ak.flatten(isolation)

ratio1=data2/data1

#Get the number of events with muons
not_empty_count1 = ak.sum(ak.num(pt_events) > 0)

nmin2=0
nmax2=5000
pt_events=MuonTree_ZeroBias["LVL1Muon_et"].array()[nmin2:nmax2]
eta_events=MuonTree_ZeroBias["LVL1Muon_eta"].array()[nmin2:nmax2]
phi_events=MuonTree_ZeroBias["LVL1Muon_phi"].array()[nmin2:nmax2]

data1=ak.flatten(pt_events)*1000 #multiply by 1000 to convert to MeV
isolation =muon_isolation_all_events(MuonTree_Zmumu, eta_events, phi_events, dr_min, dr_max, [nmin2,nmax2], 500)
data2=ak.flatten(isolation)

ratio2=data2/data1

#Get the number of events with muons
not_empty_count2 = ak.sum(ak.num(pt_events) > 0)

# %%
coolplot([ratio1, ratio2], np.linspace(-0.2,6,20),["#0072B2","#FD0000"],[fr"Z $\longrightarrow \mu \mu$ (LVL1 data),events={not_empty_count1}" ,fr"Zero Bias (LVL1 data), events={not_empty_count2}"],
        "Isolation / Transverse energy","Counts",
        fr"$\left(\frac{{E_{{t,\;iso}}}}{{E_t}}\right)$ ratio histogram, $\Delta R$ = [{dr_min}, {dr_max}]")


