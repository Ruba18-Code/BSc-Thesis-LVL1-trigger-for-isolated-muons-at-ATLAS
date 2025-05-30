
from my_functions import *

file= uproot.open("/home/ruben/Escritorio/BachelorThesisRuben/Data/Muon_trees.root") #opening the Root file with Uproot 
MuonTree_ZeroBias=file["MuonTree_ZeroBias;1"]
MuonTree_Zmumu=file["MuonTree_Zmumu;1"]


"""
Plots I want to make: 
1. Z mumu vs zero bias
2. All quality vs quality 0 for z mumu
3. Z mumu vs zero bias 0 quality
4. Different dr upper limits and lower limits

5.Repeat for LVL1 and offline
"""

"""
OFFLINE MUON PLOTS
"""

"""
Z -> mu mu vs Zero Bias Offline 
"""

#Assign eta and phi variables Zero Bias
muon_eta_all=MuonTree_ZeroBias["muon_eta"].array()
muon_phi_all=MuonTree_ZeroBias["muon_phi"].array()
nmin1=0
nmax1=2000

#Check how many events are not empty
non_empty_count1 = ak.sum(ak.num(muon_eta_all[nmin1:nmax1]) > 0)

#Compute the isolation and prepare it for plotting 
res=muon_isolation_all_events(MuonTree_ZeroBias,muon_eta_all,muon_phi_all,0.0,0.4,[nmin1,nmax1],500)
data1=ak.flatten(res)

#Assign eta and phi variables Z mu mu
muon_eta_all=MuonTree_Zmumu["muon_eta"].array()
muon_phi_all=MuonTree_Zmumu["muon_phi"].array()
nmin2=0
nmax2=300
#Check how many events are not empty
non_empty_count2 = ak.sum(ak.num(muon_eta_all[nmin2:nmax2]) > 0)

#Compute the isolation and prepare it for plotting 
res=muon_isolation_all_events(MuonTree_Zmumu,muon_eta_all,muon_phi_all,0.0,0.4,[nmin2,nmax2],100)
data2=ak.flatten(res)

#Perform the plot
colors=['b','r']
labels=[fr"Zero Bias offline muons, events={non_empty_count1}",
        fr"Z $\longrightarrow \mu \mu$ offline muons, events={non_empty_count2}"]

coolplot([data1,data2],np.linspace(0,30000,25),colors,labels,"Transverse energy (MeV)","Counts",
         r"Energy histogram isolated muons - Offline data - $\Delta R$ threshold $=$ [0,0.4]")

print("The mean values are:\n 1.Zero Bias:", int(ak.mean(data1)), "MeV \n" r" 2.Z->μμ:", int(ak.mean(data2)), "MeV\n")

"""
Z -> mu mu vs ZeroBias offline QUALITY 0
"""

#Assign eta and phi variables Zero Bias
muon_eta_all=MuonTree_ZeroBias["muon_eta"].array()
muon_phi_all=MuonTree_ZeroBias["muon_phi"].array()
nmin1=0
nmax1=2000

#Check how many events are not empty
non_empty_count1 = ak.sum(ak.num(muon_eta_all[nmin1:nmax1]) > 0)

#Compute the isolation and prepare it for plotting 
res=muon_isolation_all_events(MuonTree_ZeroBias,muon_eta_all,muon_phi_all,0.0,0.4,[nmin1,nmax1],500)
data1=ak.flatten(res)

#Assign eta and phi variables Z mu mu
muon_eta_all=quality_selector(MuonTree_Zmumu["muon_quality"].array(),MuonTree_Zmumu["muon_eta"].array(),0)
muon_phi_all=quality_selector(MuonTree_Zmumu["muon_quality"].array(),MuonTree_Zmumu["muon_phi"].array(),0)
nmin2=0
nmax2=300
#Check how many events are not empty
non_empty_count2 = ak.sum(ak.num(muon_eta_all[nmin2:nmax2]) > 0)

#Compute the isolation and prepare it for plotting 
res=muon_isolation_all_events(MuonTree_Zmumu,muon_eta_all,muon_phi_all,0.0,0.4,[nmin2,nmax2],100)
data2=ak.flatten(res)

#Perform the plot
colors=['b','r']
labels=[fr"Zero Bias offline muons, events={non_empty_count1}",
        fr"Z $\longrightarrow \mu \mu$ offline muons (QUALITY=0), events={non_empty_count2}"]

coolplot([data1,data2],np.linspace(0,30000,25),colors,labels,"Transverse energy (MeV)","Counts",
         r"Energy histogram isolated muons - Offline data - $\Delta R$ threshold $=$ [0,0.4]"
         r" - Z $\longrightarrow \mu \mu$ quality = 0")

print("The mean values are:\n 1.Zero Bias:", int(ak.mean(data1)), "MeV \n" r" 2.Z->μμ zero quality:", int(ak.mean(data2)), "MeV\n")

"""
All quality vs quality 0 for Z ->  mumu
"""

#Assign eta and phi variables Zmumu
muon_eta_all=MuonTree_Zmumu["muon_eta"].array()
muon_phi_all=MuonTree_Zmumu["muon_phi"].array()
nmin1=0
nmax1=300

#Check how many events are not empty
non_empty_count1 = ak.sum(ak.num(muon_eta_all[nmin1:nmax1]) > 0)

#Compute the isolation and prepare it for plotting 
res=muon_isolation_all_events(MuonTree_ZeroBias,muon_eta_all,muon_phi_all,0.0,0.4,[nmin1,nmax1],100)
data1=ak.flatten(res)

#Assign eta and phi variables Z mu mu
muon_eta_all=quality_selector(MuonTree_Zmumu["muon_quality"].array(),MuonTree_Zmumu["muon_eta"].array(),0)
muon_phi_all=quality_selector(MuonTree_Zmumu["muon_quality"].array(),MuonTree_Zmumu["muon_phi"].array(),0)
nmin2=0
nmax2=300
#Check how many events are not empty
non_empty_count2 = ak.sum(ak.num(muon_eta_all[nmin2:nmax2]) > 0)

#Compute the isolation and prepare it for plotting 
res=muon_isolation_all_events(MuonTree_Zmumu,muon_eta_all,muon_phi_all,0.0,0.4,[nmin2,nmax2],100)
data2=ak.flatten(res)

#Perform the plot
colors=['b','r']
labels=[fr"Z $\longrightarrow \mu \mu$ offline muons (ANY QUALITY), events={non_empty_count1}",
        fr"Z $\longrightarrow \mu \mu$ offline muons (QUALITY=0), events={non_empty_count2}"]

coolplot([data1,data2],np.linspace(0,30000,25),colors,labels,"Transverse energy (MeV)","Counts",
         r"Energy histogram isolated muons - Offline data - $\Delta R$ threshold $=$ [0,0.4]")

print("The mean values are:\n 1.Z->μμ any quality:", int(ak.mean(data1)), "MeV \n" r" 2.Z->μμ zero quality"
      , int(ak.mean(data2)), "MeV\n")

"""""""""""""""""""""""""""""""""
DIFFERENT DELTA R LIMITS - CHANGING LOWER LIMITS
"""""""""""""""""""""""""""""""""

#Assign eta and phi variables Zmumu
muon_eta_all=quality_selector(MuonTree_Zmumu["muon_quality"].array(),MuonTree_Zmumu["muon_eta"].array(),0)
muon_phi_all=quality_selector(MuonTree_Zmumu["muon_quality"].array(),MuonTree_Zmumu["muon_phi"].array(),0)
nmin1=0
nmax1=500

#Check how many events are not empty
non_empty_count1 = ak.sum(ak.num(muon_eta_all[nmin1:nmax1]) > 0)


#Compute the isolation and prepare it for plotting 

lower_dr1=0
upper_dr1=0.4
res=muon_isolation_all_events(MuonTree_ZeroBias,muon_eta_all,muon_phi_all,lower_dr1,upper_dr1,[nmin1,nmax1],100)
data1=ak.flatten(res)

lower_dr2=0.1
upper_dr2=0.4
res=muon_isolation_all_events(MuonTree_ZeroBias,muon_eta_all,muon_phi_all,lower_dr2,upper_dr2,[nmin1,nmax1],100)
data2=ak.flatten(res)

lower_dr3=0.2
upper_dr3=0.4
res=muon_isolation_all_events(MuonTree_ZeroBias,muon_eta_all,muon_phi_all,lower_dr3,upper_dr3,[nmin1,nmax1],100)
data3=ak.flatten(res)

lower_dr4=0.35
upper_dr4=0.4
res=muon_isolation_all_events(MuonTree_ZeroBias,muon_eta_all,muon_phi_all,lower_dr4,upper_dr4,[nmin1,nmax1],100)
data4=ak.flatten(res)

#Perform the plot
colors = ["#0072B2", "#FD0000", "#FFEE00", "#00FFBB"]
labels=[fr"Z $\longrightarrow \mu \mu$ (QUALITY=0), $\Delta R$ threshold = [{lower_dr1},{upper_dr1}],"
        fr"events={non_empty_count1}",
        fr"Z $\longrightarrow \mu \mu$ (QUALITY=0), $\Delta R$ threshold = [{lower_dr2},{upper_dr2}],"
        fr"events={non_empty_count1}",
        fr"Z $\longrightarrow \mu \mu$ (QUALITY=0), $\Delta R$ threshold = [{lower_dr3},{upper_dr3}],"
        fr"events={non_empty_count1}",
        fr"Z $\longrightarrow \mu \mu$ (QUALITY=0), $\Delta R$ threshold = [{lower_dr4},{upper_dr4}],"
        fr"events={non_empty_count1}"]

coolplot([data1,data2,data3,data4],np.linspace(0,30000,20),colors,labels,"Transverse energy (MeV)","Counts",
         r"Energy histogram isolated muons - Offline data - Comparing different $\Delta R$")


"""""""""""""""""""""""""""""""""
DIFFERENT DELTA R LIMITS - CHANGING UPPER LIMITS
"""""""""""""""""""""""""""""""""

#Assign eta and phi variables Zmumu
muon_eta_all=quality_selector(MuonTree_Zmumu["muon_quality"].array(),MuonTree_Zmumu["muon_eta"].array(),0)
muon_phi_all=quality_selector(MuonTree_Zmumu["muon_quality"].array(),MuonTree_Zmumu["muon_phi"].array(),0)
nmin1=0
nmax1=500

#Check how many events are not empty
non_empty_count1 = ak.sum(ak.num(muon_eta_all[nmin1:nmax1]) > 0)


#Compute the isolation and prepare it for plotting 

lower_dr1=0
upper_dr1=0.15
res=muon_isolation_all_events(MuonTree_ZeroBias,muon_eta_all,muon_phi_all,lower_dr1,upper_dr1,[nmin1,nmax1],100)
data1=ak.flatten(res)

lower_dr2=0
upper_dr2=0.4
res=muon_isolation_all_events(MuonTree_ZeroBias,muon_eta_all,muon_phi_all,lower_dr2,upper_dr2,[nmin1,nmax1],100)
data2=ak.flatten(res)

lower_dr3=0
upper_dr3=0.8
res=muon_isolation_all_events(MuonTree_ZeroBias,muon_eta_all,muon_phi_all,lower_dr3,upper_dr3,[nmin1,nmax1],100)
data3=ak.flatten(res)

lower_dr4=0
upper_dr4=np.sqrt(2)
res=muon_isolation_all_events(MuonTree_ZeroBias,muon_eta_all,muon_phi_all,lower_dr4,upper_dr4,[nmin1,nmax1],100)
data4=ak.flatten(res)

#Perform the plot
colors = ["#0072B2", "#FD0000", "#FFEE00", "#00FFBB"]
labels=[fr"Z $\longrightarrow \mu \mu$ (QUALITY=0), $\Delta R$ threshold = [{lower_dr1},{upper_dr1}],"
        fr"events={non_empty_count1}",
        fr"Z $\longrightarrow \mu \mu$ (QUALITY=0), $\Delta R$ threshold = [{lower_dr2},{upper_dr2}],"
        fr"events={non_empty_count1}",
        fr"Z $\longrightarrow \mu \mu$ (QUALITY=0), $\Delta R$ threshold = [{lower_dr3},{upper_dr3}],"
        fr"events={non_empty_count1}",
        fr"Z $\longrightarrow \mu \mu$ (QUALITY=0), $\Delta R$ threshold = [{lower_dr4},{upper_dr4}],"
        fr"events={non_empty_count1}"]

coolplot([data1,data2,data3,data4],np.linspace(0,50000,25),colors,labels,"Transverse energy (MeV)","Counts",
         r"Energy histogram isolated muons - Offline data - Comparing different $\Delta R$")

"""
LVL1 MUON PLOTS
"""

"""
Z -> mu mu vs Zero Bias LVL1 
"""

#Assign eta and phi variables Zero Bias
muon_eta_all=MuonTree_ZeroBias["LVL1Muon_eta"].array()
muon_phi_all=MuonTree_ZeroBias["LVL1Muon_phi"].array()
nmin1=0
nmax1=5000

#Check how many events are not empty
non_empty_count1 = ak.sum(ak.num(muon_eta_all[nmin1:nmax1]) > 0)

#Compute the isolation and prepare it for plotting 
res=muon_isolation_all_events(MuonTree_ZeroBias,muon_eta_all,muon_phi_all,0.0,0.4,[nmin1,nmax1],500)
data1=ak.flatten(res)

#Assign eta and phi variables Z mu mu
muon_eta_all=MuonTree_Zmumu["LVL1Muon_eta"].array()
muon_phi_all=MuonTree_Zmumu["LVL1Muon_phi"].array()
nmin2=0
nmax2=300
#Check how many events are not empty
non_empty_count2 = ak.sum(ak.num(muon_eta_all[nmin2:nmax2]) > 0)

#Compute the isolation and prepare it for plotting 
res=muon_isolation_all_events(MuonTree_Zmumu,muon_eta_all,muon_phi_all,0.0,0.4,[nmin2,nmax2],100)
data2=ak.flatten(res)

#Perform the plot
colors=['b','r']
labels=[fr"Zero Bias LVL1 muons, events={non_empty_count1}",
        fr"Z $\longrightarrow \mu \mu$ LVL1 muons, events={non_empty_count2}"]

coolplot([data1,data2],np.linspace(0,30000,25),colors,labels,"Transverse energy (MeV)","Counts",
         r"Energy histogram isolated muons - LVL1 data - $\Delta R$ threshold $=$ [0,0.4]")

print("The mean values are:\n 1.Zero Bias:", int(ak.mean(data1)), "MeV \n" r" 2.Z->μμ:", int(ak.mean(data2)), "MeV\n")


"""""""""""""""""""""""""""""""""
DIFFERENT DELTA R LIMITS - CHANGING LOWER LIMITS
"""""""""""""""""""""""""""""""""

#Assign eta and phi variables Zmumu
muon_eta_all=MuonTree_Zmumu["LVL1Muon_eta"].array()
muon_phi_all=MuonTree_Zmumu["LVL1Muon_phi"].array()
nmin1=0
nmax1=500

#Check how many events are not empty
non_empty_count1 = ak.sum(ak.num(muon_eta_all[nmin1:nmax1]) > 0)


#Compute the isolation and prepare it for plotting 

lower_dr1=0
upper_dr1=0.4
res=muon_isolation_all_events(MuonTree_ZeroBias,muon_eta_all,muon_phi_all,lower_dr1,upper_dr1,[nmin1,nmax1],100)
data1=ak.flatten(res)

lower_dr2=0.1
upper_dr2=0.4
res=muon_isolation_all_events(MuonTree_ZeroBias,muon_eta_all,muon_phi_all,lower_dr2,upper_dr2,[nmin1,nmax1],100)
data2=ak.flatten(res)

lower_dr3=0.2
upper_dr3=0.4
res=muon_isolation_all_events(MuonTree_ZeroBias,muon_eta_all,muon_phi_all,lower_dr3,upper_dr3,[nmin1,nmax1],100)
data3=ak.flatten(res)

lower_dr4=0.35
upper_dr4=0.4
res=muon_isolation_all_events(MuonTree_ZeroBias,muon_eta_all,muon_phi_all,lower_dr4,upper_dr4,[nmin1,nmax1],100)
data4=ak.flatten(res)

#Perform the plot
colors = ["#0072B2", "#FD0000", "#FFEE00", "#00FFBB"]
labels=[fr"Z $\longrightarrow \mu \mu$ LVL1 muons, $\Delta R$ threshold = [{lower_dr1},{upper_dr1}],"
        fr"events={non_empty_count1}",
        fr"Z $\longrightarrow \mu \mu$ LVL1 muons, $\Delta R$ threshold = [{lower_dr2},{upper_dr2}],"
        fr"events={non_empty_count1}",
        fr"Z $\longrightarrow \mu \mu$ LVL1 muons, $\Delta R$ threshold = [{lower_dr3},{upper_dr3}],"
        fr"events={non_empty_count1}",
        fr"Z $\longrightarrow \mu \mu$ LVL1 muons, $\Delta R$ threshold = [{lower_dr4},{upper_dr4}],"
        fr"events={non_empty_count1}"]

coolplot([data1,data2,data3,data4],np.linspace(0,30000,20),colors,labels,"Transverse energy (MeV)","Counts",
         r"Energy histogram isolated muons - LVL1 data - Comparing different $\Delta R$")


"""""""""""""""""""""""""""""""""
DIFFERENT DELTA R LIMITS - CHANGING UPPER LIMITS
"""""""""""""""""""""""""""""""""

#Assign eta and phi variables Zmumu
muon_eta_all=MuonTree_Zmumu["LVL1Muon_eta"].array()
muon_phi_all=MuonTree_Zmumu["LVL1Muon_phi"].array()
nmin1=0
nmax1=500

#Check how many events are not empty
non_empty_count1 = ak.sum(ak.num(muon_eta_all[nmin1:nmax1]) > 0)


#Compute the isolation and prepare it for plotting 

lower_dr1=0
upper_dr1=0.15
res=muon_isolation_all_events(MuonTree_ZeroBias,muon_eta_all,muon_phi_all,lower_dr1,upper_dr1,[nmin1,nmax1],100)
data1=ak.flatten(res)

lower_dr2=0
upper_dr2=0.4
res=muon_isolation_all_events(MuonTree_ZeroBias,muon_eta_all,muon_phi_all,lower_dr2,upper_dr2,[nmin1,nmax1],100)
data2=ak.flatten(res)

lower_dr3=0
upper_dr3=0.8
res=muon_isolation_all_events(MuonTree_ZeroBias,muon_eta_all,muon_phi_all,lower_dr3,upper_dr3,[nmin1,nmax1],100)
data3=ak.flatten(res)

lower_dr4=0
upper_dr4=np.sqrt(2)
res=muon_isolation_all_events(MuonTree_ZeroBias,muon_eta_all,muon_phi_all,lower_dr4,upper_dr4,[nmin1,nmax1],100)
data4=ak.flatten(res)

#Perform the plot
colors = ["#0072B2", "#FD0000", "#FFEE00", "#00FFBB"]
labels=[fr"Z $\longrightarrow \mu \mu$ LVL1 muons, $\Delta R$ threshold = [{lower_dr1},{upper_dr1}],"
        fr"events={non_empty_count1}",
        fr"Z $\longrightarrow \mu \mu$ LVL1 muons, $\Delta R$ threshold = [{lower_dr2},{upper_dr2}],"
        fr"events={non_empty_count1}",
        fr"Z $\longrightarrow \mu \mu$ LVL1 muons, $\Delta R$ threshold = [{lower_dr3},{upper_dr3}],"
        fr"events={non_empty_count1}",
        fr"Z $\longrightarrow \mu \mu$ LVL1 muons, $\Delta R$ threshold = [{lower_dr4},{upper_dr4}],"
        fr"events={non_empty_count1}"]

coolplot([data1,data2,data3,data4],np.linspace(0,50000,25),colors,labels,"Transverse energy (MeV)","Counts",
         r"Energy histogram isolated muons - LVL1 data - Comparing different $\Delta R$")
