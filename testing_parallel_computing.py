# %%
from my_functions import *



# %%
file= uproot.open("/home/ruben/Escritorio/BachelorThesisRuben/Data/Muon_trees.root") #opening the Root file with Uproot 

MuonTree_ZeroBias=file["MuonTree_ZeroBias;1"]


# %%
def dr_threshold_boolean_mask_event(event_dr,threshold):

    """"
    Takes an array containing the delta r values for an event and its respective jTower,
    then it starts to iterate over al muons in the event and generates a boolean mask containing 'True'
    for all delta_r values that are below a certain threshold.

    By using NumPy we're able to do this in a vectorized way (meaning: operating over all the array at once, instead of looping), which is
    shorter and apparently more resource efficient

    The idea is to apply this boolean mask to the energy vector, that will select the energy elements that are associated with the muon
    """
    event_dr = np.array(event_dr)  # ensures it's a NumPy array
    return event_dr < threshold

# %%
def muon_isolation_one_event(muon_eta_event, muon_phi_event, jTower_eta_event, jTower_phi_event, jTower_et_event,threshold):

    """
    Inputs:
        -Array containing the values of LVL1 muon pseudorapidity for an event (muon_eta_event)
        -Array containing the values of LVL1 muon phi for an event (muon_phi_event)
        -Array containing the values of LVL1 muon transverse energy in MeV for an event (muon_et_event)

        -Array containing the values of jTower pseudorapidity for an event (jTower_eta_event)
        -Array containing the values of jTower phi for an event (jTower_phi_event)
        -Array containing the values of jTower transverse energy in MeV for an event (jTower_et_event)

        -A scalar that indicates the threshold for delta r (threshold)

    It will return an array containing the isolated muon energy of each muon inside the event.

    To do so, it created pairs of (muon_eta,jTower_eta), (muon_phi,jTower_phi) for each muon and each jTower element, used to
    compute the delta r values, and finally it creates a boolean mask that depends on the threshold and selects the valid calorimeter energies.
    """

    isolated_energy_event=[]

    for (eta, phi) in zip(muon_eta_event, muon_phi_event):
        #Generate pairs of [muon,jTower1],[muon,jTower2]...
        #That's an array where the left column is the value of a muon (always the same) and the right column contains all the jTower values
        #associated with such muon
        jTower_muon_eta_pairs=[(eta, stuff) for stuff in jTower_eta_event]
        jTower_muon_phi_pairs=[(phi, stuff) for stuff in jTower_phi_event]

        #Compute delta r with the (eta,phi) pairs
        dr_jTower_muon=delta_r( jTower_muon_eta_pairs,jTower_muon_phi_pairs)

        #Create a boolean mask that will be 'True' only if the computed dr is smaller or equal than the threshold
        mask=dr_threshold_boolean_mask_event(dr_jTower_muon,threshold)

        #Apply the mask to the energy vector. This will select only the energies corresponding to the jTower pixels marked as TRUE
        result=jTower_et_event[mask]

        #Take the negative values out (they appear due to technical features of the detector)
        #Sum up the results to get an estimate of the isolated muon energy in MeV
        T=sum(result[result >= 0])
        isolated_energy_event.append(T)

    return(isolated_energy_event)

# %%
# %%
import dask
from dask import delayed
import dask.array as da
from tqdm import tqdm

import dask
from dask import delayed
from dask.diagnostics import ProgressBar
from tqdm import tqdm

def muon_isolation_all_events(MuonTree_ZeroBias, threshold, batch_size=100):
    # Extract all the muon eta and phi values
    muon_eta_all_events = MuonTree_ZeroBias["LVL1Muon_eta"].array()
    muon_phi_all_events = MuonTree_ZeroBias["LVL1Muon_phi"].array()

    # List to store delayed computations
    res = []

    # Define a function to process each event
    def process_event(n):
        # Extract each event's data
        muon_eta_event = muon_eta_all_events[n]
        muon_phi_event = muon_phi_all_events[n]

        # Extract jTower data for the same event
        jTower_eta_event, jTower_phi_event, jTower_et_event = jTower_selector(
            MuonTree_ZeroBias, ["jTower_eta", "jTower_phi", "jTower_et_MeV"], n
        )

        # Compute the isolation for this event
        isol_event = muon_isolation_one_event(
            muon_eta_event, muon_phi_event, jTower_eta_event, jTower_phi_event, jTower_et_event, threshold
        )

        return isol_event

    # Get total number of events
    total_events = len(muon_eta_all_events)

    # Process the events in smaller batches
    for start in tqdm(range(0, 1000, batch_size)):  # Loop in steps of batch_size
        end = min(start + batch_size, total_events)  # Ensure we don't exceed total events
        batch = [delayed(process_event)(n) for n in range(start, end)]  # Create delayed tasks for this batch
        res.extend(batch)

        # Compute this batch of tasks, wait for results to free up memory
        with ProgressBar():
            dask.compute(*batch)  # Compute all delayed tasks in this batch

    return res




# %%



print(muon_isolation_all_events(MuonTree_ZeroBias,0.4))

