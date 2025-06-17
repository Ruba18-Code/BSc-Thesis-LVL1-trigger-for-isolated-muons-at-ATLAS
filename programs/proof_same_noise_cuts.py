from my_functions import *

file= uproot.open("/home/ruben/Escritorio/BachelorThesisRuben/Data/Muon_trees.root") #opening the Root file with Uproot 
MuonTree_ZeroBias=file["MuonTree_ZeroBias;1"]
MuonTree_Zmumu=file["MuonTree_Zmumu;1"]


"""
This script is used to prove that the noise cuts are the same for all events and there's no need to compute the cuts for each event.
I can compute the cuts once (for the first event) and then use the same cuts for all events. This is much faster.
"""
#Choose the range of events to plot
nmin=0
nmax=100

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

def muon_isolation_all_events2(tree,muon_eta_all,muon_phi_all, lower_threshold, upper_threshold, event_range, batch_size=10000, get_mask=False):
    """
    This function computes muon isolation for all events in batches of a certain size to avoid crashing the computer.
    
    Inputs:
        tree: select a tree (MuonTree_ZeroBias or MuonTree_Zmumu)
        muon_eta/phi_all: awkward array containing the data of a Muon Tree (example: muon_eta_all=MuonTree_ZeroBias["LVLMuon_eta"].array())
        threshold: float that sets the maximum acceptable value for delta R (typically 0.4)
        event_range: tuple of [start_event, end_event]. To use the full dataset write [0,len(muon_eta_all)]
        batch_size: int that sets the number of events to process at once (shouldn't be larger than 5.000-10.000 for a PC with 8gb of RAM)
        get_mask: boolean that indicates if the mask should be returned or not, False by default
    Returns:
        List of isolated muon energies per event
    """
    #First compute the cuts for the jTower energies
    #jTower_cuts=jTower_assign_cuts(tree)

    start_event, end_event = event_range
    res = []
    masks=[]

    # Process in batches to avoid memory overload (crash due to lack of RAM)

    total_events = end_event - start_event

    #tqdm is used to get a progress bar that estimates the remaining time 
    #batch_stop is written like that to prevent going out of range
    for batch_start in tqdm(range(0, total_events, batch_size), desc="muon_isolation_all_events: Computing muon isolation", leave=False):
        batch_stop = min(batch_start + batch_size, total_events)

        # Load jTower data only for the current batch
        jTower = tree.arrays(
            ["jTower_eta", "jTower_phi", "jTower_et_MeV"],
            entry_start=start_event + batch_start,
            entry_stop=start_event + batch_stop
        )
        #assign eta, phi and et
        jTower_eta_batch = jTower["jTower_eta"]
        jTower_phi_batch = jTower["jTower_phi"]
        jTower_et_batch  = jTower["jTower_et_MeV"]

        #get the muon data for the current batch
        muon_eta_batch = muon_eta_all[batch_start:batch_stop]
        muon_phi_batch = muon_phi_all[batch_start:batch_stop]

        # Loop over events in a batch
        for i in range(len(muon_eta_batch)):
            event_index = start_event + batch_start + i
            jTower_cuts=ak.flatten(jTower_assign_cuts(tree,event_index,event_index+1))
            muon_eta_event = muon_eta_batch[i]
            muon_phi_event = muon_phi_batch[i]

            #Check if the event is empty before starting the calculations. This can reduce the computing time a lot if
            #our data includes many empty events 
            if len(muon_eta_event) == 0:
                res.append([])
            else:
                #If the event is not empty compute the isolation
                jTower_eta_event = jTower_eta_batch[i]
                jTower_phi_event = jTower_phi_batch[i]
                jTower_et_event  = jTower_et_batch[i]

                isol_event, mask= muon_isolation_one_event(
                    muon_eta_event, muon_phi_event,
                    jTower_eta_event, jTower_phi_event, jTower_et_event,
                    lower_threshold,upper_threshold,jTower_cuts,get_mask)
                #and append it to the result
                res.append(isol_event)
                #if get_mask is True, append the mask to the result
                if get_mask==True:
                    masks.append(mask)
    #if get_mask is True, return the result and the mask, otherwise just the result
    if get_mask==True:
        return res, masks
    else:
        return res
    
res1=muon_isolation_all_events(MuonTree_Zmumu,Zmumu_eta,Zmumu_phi, 0.0, 0.4, [0,nmax])
res2=muon_isolation_all_events2(MuonTree_Zmumu,Zmumu_eta,Zmumu_phi, 0.0, 0.4, [0,nmax])

print(res1==res2)