# %%
import uproot
import numpy as np #math and science package
import scipy as sp #math and science package
import awkward as ak #root files are usually awkward arrays 
import matplotlib.pyplot as plt #plot stuff
import itertools
from tqdm import tqdm #generated progress bars and estimates the computation time
from numba import njit #numba helps to speed up calculations if the code is written in a certain way

def histogram(data, nbins, x_range,  x_label, y_label, title_label):
    plt.hist(data, bins=nbins, range=x_range, color=plt.cm.viridis(0.9), edgecolor= 'black', density=True, alpha=0.7)   #plot the histogram of the data, for example 100 bins between -3 and 3
    plt.xlabel(x_label) #to wrtie LaTeX letters, use "r" at the beginning 
    plt.ylabel(y_label)
    plt.title(title_label)
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.show()

def histogram2(data1, data2, nbins, x_range,  x_label, y_label, title_label, label1, label2):

    #This function aims to plot two histograms at once. We need to provide two datasets together with some parameters like number of 
    #bins and label

    if not isinstance(data1, list):
        data1=ak.to_numpy(ak.flatten(data1))
    if not isinstance(data2, list):
        data2=ak.to_numpy(ak.flatten(data2))
    
    plt.hist(data1, bins=nbins, range=x_range, edgecolor= 'green', density=True, label=label1, alpha=1,histtype='step') #density=True normalises
    plt.hist(data2, bins=nbins, range=x_range, edgecolor= 'blue', density=True, label=label2, alpha=1,histtype='step')
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(title_label)
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.legend()
    plt.show()

def histogram2errors(data1, data2, nbins, x_range, x_label, y_label, title_label, label1, label2):

    #This fuctions plots two histograms at once and adds the errorbars

    # Flatten awkward arrays
    data1 = ak.flatten(data1, axis=None)
    data2 = ak.flatten(data2, axis=None)

    # Compute histograms
    y1, binEdges = np.histogram(data1, bins=nbins, range=x_range, density=True)
    y2, _ = np.histogram(data2, bins=nbins, range=x_range, density=True)
    
    # Bin centers and width
    bincenters = 0.5 * (binEdges[1:] + binEdges[:-1])
    bin_width = binEdges[1] - binEdges[0]
    width = bin_width * 0.4  # Width for each bar (40% for each group, leaving 20% gap)

    # Calculate statistical uncertainties (Poisson errors assuming density=True)
    counts1, _ = np.histogram(data1, bins=nbins, range=x_range)
    counts2, _ = np.histogram(data2, bins=nbins, range=x_range)
    norm1 = np.sum(counts1)
    norm2 = np.sum(counts2)
    yerr1 = np.sqrt(counts1) / norm1 / bin_width
    yerr2 = np.sqrt(counts2) / norm2 / bin_width

    # Plot with error bars
    plt.bar(bincenters - width/2, y1, width=width, edgecolor='black', color='r', alpha=0.8, label=label1, yerr=yerr1, capsize=2)
    plt.bar(bincenters + width/2, y2, width=width,edgecolor='black',color='b', alpha=0.8, label=label2, yerr=yerr2, capsize=2)

    # Labels and legend
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(title_label)
    plt.legend()
    plt.xlim(x_range)
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.show()

def coolplot(data,bins,colors,labels,x_label,y_label,title):

    #This function is designed to create a comparative histogram of N sets of data with their respective errorbars included.
    #Introduce a data vector, that contains one column for each dataset that has to be plotted data=[[set_a],[set_b],...]
    #bins should be a linspace containing np.linspace(start:end, number_of_bins)
    #colors is a vector of colors
    #labels is a vector of labels for the legend
    #x and y_label are NOT vectors: name of x and y axis
    #same for title 

    #All outsider values will appear concentrated on the first and last bins

    bins_overflow = np.concatenate([[-np.inf],bins[1:-1],[np.inf]]) 

    #Create empty array

    hists=[]
    error=[]
    hists_norm=[]
    error_norm=[]

    #Compute the histogram and error for every dataset
    for i in range(len(data)):
        hists.append(np.histogram(ak.flatten(data[i], axis=None),bins=bins_overflow)[0])
        error.append(np.sqrt(hists[i]))

    #Plot according to some preferences

    for i in range(len(data)):
        error[i]=error[i]/np.sum(hists[i])
        hists[i]=hists[i]/np.sum(hists[i])

        #Plot function
        plt.step(
            bins[:-1],
            hists[i],
            where='mid',  # or 'post', or 'pre' depending on visual preference
            color=colors[i],
            label=labels[i],
            alpha=0.5)

        #Plot errorbars
        plt.errorbar(x=bins[:-1],
            y=hists[i],
            yerr=error[i],
            color='black',
            fmt='o',
            markersize=1)
    
    #cool parameters + plot the title and labels

    plt.tight_layout()
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(title)
    plt.show()

#--------------------------------------------------------------------------------

def ak_element_lenght_counter(data):

    #This function takes an awkward array and returns an array that contains the lenght of every event
    #Example: data=([1,1],[2,2,2,2]) will return (2,4)
    return(ak.num(data))

def plot_number_elements_per_event(data,bins,colors,labels,x_label,y_label,title):
    
    #The following function combines ak_element_length_counter and coolplot. It plots an histogram with the number of elements per event

    event=[]
    
    for i in range(len(data)):
        a=ak_element_lenght_counter(data[i])
        event.append(a)

    #-------PLOT-----------------------------
    #            
    coolplot(event,bins,colors,labels,x_label,y_label,title)

#-----------------------------------------------------------------------------------------------

def invariant_mass_pair(pt1, eta1, phi1, pt2, eta2, phi2):
    """
    Compute the invariant mass of a particle pair system given their pt, eta, and phi.
    
    Inputs:
        -pt1, eta1, phi1: floats, transverse momentum, pseudorapidity, and azimuthal angle of particle 1
        -pt2, eta2, phi2: floats, transverse momentum, pseudorapidity, and azimuthal angle of particle 2

    Returns:
        -float: Invariant mass of the system in MeV
    """

    # Compute cartesian momenta
    px1, py1 = pt1 * np.cos(phi1), pt1 * np.sin(phi1)
    pz1 = pt1 * np.sinh(eta1)
    E1  = pt1 * np.cosh(eta1)

    px2, py2 = pt2 * np.cos(phi2), pt2 * np.sin(phi2)
    pz2 = pt2 * np.sinh(eta2)
    E2  = pt2 * np.cosh(eta2)

    # Total momentum and energy
    px, py, pz = px1 + px2, py1 + py2, pz1 + pz2
    E = E1 + E2

    # Invariant mass squared
    mass_squared = E**2 - (px**2 + py**2 + pz**2)

    #Check if it's positive for safety and return
    return np.sqrt(mass_squared) if mass_squared > 0 else 0.0


def invariant_mass_pair_selector(pt,eta,phi):

    """
    This function is designed to take pt (in MeV!), eta and phi from a Z-> mu mu event. If the event involves more than 2 muons, that means
    that there are 'fake' muons included in our event. The idea is to compute the invariant mass for all possible pairs
    and select which is the 'True' pair by comparint the computed invariant mass with the Z boson invariant mass. The pair
    that gets closer will be considered the best candidate

    Inputs:
        -pt, eta and phi of the event
    
    Returns:
        -pt, eta and phi of the chosen pair
        -the invariant mass of the best pair
    """

    #Let's hardcode the theoretical value of the Z boson invariant mass
    ztheo=91188 #MeV

    #If the length of the data is not consistent return empty and None
    if not (len(pt) == len(eta) == len(phi)):
        pair=ak.Array([None,None,None])
        inv_mass=None

        return(inv_mass, pair)
    
    #If the length is consistent start selecting
    else:
        #If the event has 1 or 0 muons there's no possible pair to be formed: return empty and None
        if len(pt)<2:
            pair=ak.Array([None,None,None])
            inv_mass=None

            return(inv_mass, pair)
        
        #If the event contains 2 muons (1 pair)
        if len(pt)==2:
            #The selected pair will be of course the only pair
            pair=([pt,eta,phi])
            #And compute the invariant mass of the pair
            inv_mass=invariant_mass_pair(pt[0],eta[0],phi[0],pt[1],eta[1],phi[1])

            return(inv_mass,pair)
        #If the event contains more than 2 muons, take into account all the possible combinations
        if len(pt)>2:
            mu=[]
            #for each muon in the event, take its pt,eta,phi
            for (muon_pt,muon_eta,muon_phi) in zip(pt,eta,phi):
                #append it to a list containing pt,eta,phi for all muons in the event ([mu1],[mu2],[mu3])
                mu.append([muon_pt,muon_eta,muon_phi])
            #Using itertools create all possible pairs of muons
            pairs=list(itertools.combinations(mu,2))

            #for every muon pair
            res=[]
            diff=[]
            for mu1, mu2 in pairs:
                #unpack the variables
                pt1, eta1, phi1 = mu1
                pt2, eta2, phi2 = mu2

                #Compute the invariant mass of the pair
                inv_mass=invariant_mass_pair(pt1, eta1, phi1, pt2, eta2, phi2)
                res.append(inv_mass)
                #Compute the difference with respect to the theoretical Z mass
                pair_diff=abs(ztheo-inv_mass)
                diff.append(pair_diff)

            #Now select the best pair (meaning: pair with the smallest differencie with respect to the theoretical value)
            diff = np.array(diff)
            res = np.array(res)
            pairs = np.array(pairs)

            # Boolean mask where the diff is equal to the minimum
            mask = diff == np.min(diff)

            # Apply mask to res and pairs ([0] selects the first in the unlikely case of having two muons with the same diff)
            inv_mass= res[mask][0]
            pair = pairs[mask][0]

            return(inv_mass,pair)
        
        #else for safety
        else:
            pair=ak.Array([None,None,None])
            inv_mass=None

            return(inv_mass, pair)


def invariant_mass_all_muons(pt,eta,phi):

    """
    Introduce the arrays containing pt,eta and phi data for an event of arrays of events. The function will compute the 
    invariant mass for each event. If the event involves more than two muons, it will select the muon pair that gets
    closer to the Z boson invariant mass, assuming that this is the 'True' pair of muons between all the other candidates.

    Inputs:
        -pt, eta, phi: awkward arrays, data containing pt, eta and phi for one or multiple events involving any number of muons
    
    Returns:
        -awkward array containing the invariant mass for each event. If there's more than 2 muons, it selects the one that's closer
        to the Z boson invariant mass.

        -tqdm generates a progress bar to help preserve the user's mental health
    """

    invariant_masses=[]
    
    #Scan all events
    for i in tqdm(range(len(eta))):

        #If the event involves two muons
        if len(eta[i])==2 and len(phi[i])==2 and len(pt[i]) == 2:

            #Unpack
            eta1, eta2= eta[i]
            phi1, phi2 = phi[i]
            pt1, pt2= pt[i]

            #Compute invariant mass
            invariant_masses.append(invariant_mass_pair(pt1, eta1, phi1, pt2, eta2, phi2))

        #If it involves more than 2 muons (also check if they're the same length for safety)
        elif len(pt[i]) > 2 and len(pt[i]) == len(eta[i]) == len(phi[i]):
            m = invariant_mass_pair_selector(pt[i], eta[i], phi[i])
            if m is not None:  # Ensure that m is not None
                invariant_masses.append(m[0])
        else:
            invariant_masses.append([])

    return ak.Array(invariant_masses)

def get_all_Z_peak_pairs(pt, eta, phi):

    """This function returns the selected pairs during the Z-peak reconstruction. Meaning: if we introduce Z->mu mu events involving more
    than two muons, it will return events involving only the 2 muons, because it filters out the 'Fake' muons by computing the invariant
    mass of all possible muon pairs in an events and comparing it to the theoretical Z boson mass (all that is done using invariant_mass
    _pair_selector). The role of this functions is structuring the output of invariant_mass_pair_selector in a useful way.

    Inputs:
        -pt, eta, phi: awkward arrays, transverse momentum, pseudorapidity and azimuthal angle of the events

    Returns:
        -pairs: nested awkward array with shape [event1,event2,...] where event1=[mu1,mu2] where mu1=[pt1,eta1,phi1].
    
    """
    pairs = []
    for i in tqdm(range(len(eta))):
        #Execute invariant_mass_pair_selector for the i-th event and get the selected pair (adding [1] at the end)
        selected_pair = (invariant_mass_pair_selector(pt[i], eta[i], phi[i]))[1]

        #Check if it's None for safety
        if selected_pair is not None:
            # Check if selected_pair is structured as expected
            if len(selected_pair) == 3:
                pt_pair, eta_pair, phi_pair = selected_pair
                
                # Make sure none of these are None or empty
                if (pt_pair is not None and eta_pair is not None and phi_pair is not None):

                    #Re-structure the data in a more useful way
                    combined_pair = [list(x) for x in zip(pt_pair, eta_pair, phi_pair)]
                    pairs.append(combined_pair)
                    continue  # continue to next event
        
        # If we got here, something was missing or None - append empty list
        pairs.append([[],[]])
    
    return ak.Array(pairs)


#-----------------------------------------------------------------------------------------------

def quality_selector(quality_data, muon_data, value):

    """
    This function is designed to select only the events that match a certain quality.

    Inputs:
        -quality_data: awkward array, contains the quality data (in my case: MuonTrees_Zmumu["muon_quality"])
        -muon_data: awkward array, contains the data that we want to filter if it matches the specified quality 
        (example: MuonTrees_Zmumu["muon_eta"])
        -value: int, desired quality (0 represents the highest quality)

    Returns:
        -A filtered version of muon_data containing only events where all muons have quality==value
        -Prints a message that informs about how much data has been selected.

    """ 
    # Create a mask that is True if ALL elements in the event match the value
    event_mask = ak.all(quality_data == value, axis=1)

    # Apply the mask to select only the events where all muons have the desired quality
    selected_data = muon_data[event_mask]

    # Compute the selection percentage
    print("Only", (len(selected_data) / len(muon_data)) * 100, r"% of the data has been selected")

    return selected_data


#-----------------------------------------------------------------------------------------------
def delta_phi(phipairs):

    """
    This funtion takes an array containing pais of phi: phipairs=([phi1,phi2],[phi3,phi4],...)
    and computes the delta phi for each pair in a vectorized way. 

    It also works for just one pair.

    The calculation on the 'return' line ensures that the output is between -Pi and +Pi
    """
    #Take first column - second column
    phipairs= np.atleast_2d(phipairs)
    dphi = phipairs[:, 0] - phipairs[:, 1]

    #Return value between -Pi and +Pi
    return (dphi + np.pi) % (2 * np.pi) - np.pi

def delta_eta(etapairs):

    """
    This funtion takes an array containing pais of eta: etapairs=([eta1,eta2],[eta3,eta4],...)
    and computes the delta eta for each pair in a vectorized way. 

    It also works for just one pair.
    """

    #This line ensures that it'll work for 1D arrays [eta1,eta2] adding an extra dimension -> [[eta1,eta2]]
    etapairs = np.atleast_2d(etapairs)
    return(np.abs(etapairs[:,0] - etapairs[:,1]))

def delta_r(etapairs, phipairs):

    """
    Computes delta r (Δr) between eta and phi pairs in a vectorized way.
    This funtion takes an array containing pais of eta: etapairs=([eta1,eta2],[eta3,eta4],...)
    and also an array containing pais of phi: phipairs=([phi1,phi2],[phi3,phi4],...)

    It also works for just one pair.
    """
    return np.hypot(delta_eta(etapairs),delta_phi(phipairs))
#---------------------------------------------------------------------------

def rate_calculator(data,cut,scaling_factor):

    #This function aims to compute the rate of events that are above a certain cut
    #Example: our data is muon energy, we say cut=10000MeV, then the function will return
    #the fraction of events that are above this cut (rate) and an array with these events (aux)

    #moreover, 'scaling' allows us to scale the rate (by default between 0 and 1), choose scaling=1 if no scaling is needed
    #in my case, scaling=40*10⁶*(2340/3564)Hz since this is the frequency of filled bunch crossings at ATLAS

    #Creates a boolean mask that contains 'True' if the data is above the cut
    data_above_cut = ak.any(data >= cut, axis=-1)

    #Apply boolean mask to extract the elements related to 'True'
    aux = data[data_above_cut]
    
    # Calculate the rate
    rate = len(aux) / len(data)
    rate=rate*scaling_factor

    return rate, aux

#---------------------------------------------------------------------------------------

def jTower_handler(tree,name,xmin,xmax,binsize,chunk_size,xlabel,ylabel,title,showplot):

    #This function aims to manipulate the jTower data
    #Since these datasets are very large, it divides them into chunks (recommended value: chunk_size= 10000)
    #It also returns how many chunks have been made ('steps')
    #Finally, it calculates a histogram with the sum of all the chunk histograms and plots it if showplot is True

    #'name' must be a string of characters like: "muon_et" for instance

 #Initialize these variables
 steps=0
 hists=[]

 #Create the x axis 
 bins=np.arange(xmin,xmax,binsize)

 #Slice the dataset in chunks of size 'chunk_size' and count the amount of chunks 
 for data in tree.iterate([name], step_size=chunk_size):
       data = data[name]
       steps=steps+1
        #Compute histogram for all chunks
       hists.append(np.histogram(ak.flatten(data, axis=None),bins)[0])

 #Add the histograms together
 combined_histogram=sum(hists)

 #Plot the histogram
 if showplot == True:
      bins=np.arange(xmin,xmax-binsize,binsize) #There is one bin less, this is necessary
      plt.plot(bins,combined_histogram,color='b')
      plt.xlabel(xlabel)
      plt.ylabel(ylabel)
      plt.title(title)
      plt.grid(True,linestyle='--', alpha=0.5)
      plt.show()
      
 return(combined_histogram, steps)

def jTower_selector(tree, names, n):

    """
    Selects specific branches for a given event index.

    Inputs:
    - tree: awkward array tree (e.g. Muon_ZeroBias or Muon_Zmumu)
    - names: A string or list of strings representing the branch names (e.g. ["jTower_eta", "jTower_et_MeV"])
    - n: index of the event

    Returns:
    - A single value if 'names' is a string (just 1 name), or a list of values if 'names' is a list of strings (more than 1 name)
    """

    if isinstance(names, str): #if the input is a string (just 1 branch name has been introduced)
        return tree.arrays([names], entry_start=n, entry_stop=n + 1)[names][0]
    else: #If the input is a list of strings (more than 1 branch name introduced)
        result = tree.arrays(names, entry_start=n, entry_stop=n + 1)
        return [result[name][0] for name in names]


# %%
def jTower_muon_pair_maker(tree, name, n , val):

    """
    This function takes the jTower values of a certain tree (MuonTree_ZeroBias or MuonTree_Zmumu)
    and a certain variable 'name' ("jTower_eta","jTower_phi","jTower_et_MeV"... IT HAS TO BE A STRING " ") for a specific event denoted with the index 'n' (n-th event).
    And a scalar 'val' that should correspond to the value of eta,phi,et... of a muon in the n-th event.

    Then it creates an array containing pairs, where the left column is always 'val' and the right column contains all
    the elements inside our jTower event.

    This is useful to then insert the output inside the delta_eta and delta_phi functions, and conveniently compute all
    delta_r between the muon and the jTower.

    The 'return' line is written in that way instead of what I though initially:

    res=[]
    #Create pairs
    for stuff in jTower_array:
        tup= [val, stuff]
        res.append(tup)

    return(res)

    because it's faster (we avoid .append)

    """

    #Get the jTower array
    jTower_array=jTower_selector(tree,name,n)

    # Use list comprehension to create pairs (val, jTower_element) for each jTower_element.
    return [(val, stuff) for stuff in jTower_array]

############################################################3

def isElement_below_threshold(element, threshold):
    """
    Checks if an element (scalar or awkward array) is below a threshold.
    Returns True if any value is below the threshold.
    Returns False if the element is empty or None.
    """
    #Check if its None to avoid errors
    if element is None:
        return False
    
    #If it's not an awkward array then convert it to one
    if not isinstance(element, ak.Array):
        element = ak.Array([element])

    #Return true only if it's below threshold, False otherwise
    return ak.any(element <= threshold)


def isMuon_below_threshold(muon, threshold):
    """
    For a muon (list of elements), return a list of booleans:
    True where any element is below threshold, False otherwise.
    """
    #Check if it's none of empty
    if muon is None or len(muon) == 0:
        return ak.Array([False])
    
    #For all elements in the muon, execute the isElement function
    return ak.Array([isElement_below_threshold(el, threshold) for el in muon])

def isEvent_below_threshold(event, threshold):
    """
    For an event (list of muons), return a list of booleans:
    True where any element in a muon is below threshold.
    """

    #Check if it's None or empty
    if event is None or len(event) == 0:
        return ak.Array([False])
    
    #For all muons in an event, execute the isMuon function
    return ak.Array([isMuon_below_threshold(muon, threshold) for muon in event])

def isAll_below_threshold(data, threshold):
    """
    This function takes a dataset (awkward array, potentially nested) and a threshold (scalar)
    and returns an array with the same shape containing True or False depending on if the 
    elements of the data set are below the threshold or not.

    This has been designed to check which elements of the delta r array are below a desired distance
    """

    #If the data is None or empty
    if data is None or len(data) == 0:
        return ak.Array([False])
    
    #For all events in the data, execute the isEvent function
    return ak.Array([isEvent_below_threshold(event, threshold) for event in data])
#-----------------------------------------------------------------------------------------------
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
        T=np.sum(result[result >= 0])
        isolated_energy_event.append(T)

    return(isolated_energy_event)

def muon_isolation_all_events(tree,muon_eta_all,muon_phi_all, threshold, event_range, batch_size):
    """
    This function computes muon isolation for all events in batches of a certain size to avoid crashing the computer.
    
    Inputs:
        tree: select a tree (MuonTree_ZeroBias or MuonTree_Zmumu)
        muon_eta/phi_all: awkward array containing the data of a Muon Tree (example: muon_eta_all=MuonTree_ZeroBias["LVLMuon_eta"].array())
        threshold: float that sets the maximum acceptable value for delta R (typically 0.4)
        event_range: tuple of [start_event, end_event]. To use the full dataset write [0,len(muon_eta_all)]
        batch_size: int that sets the number of events to process at once (shouldn't be larger than 5.000-10.000 for a PC with 8gb of RAM)
    Returns:
        List of isolated muon energies per event
    """

    start_event, end_event = event_range
    res = []

    # Process in batches to avoid memory overload (crash due to lack of RAM)

    total_events = end_event - start_event

    #tqdm is used to get a progress bar that estimates the remaining time 
    #batch_stop is written like that to prevent going out of range
    for batch_start in tqdm(range(0, total_events, batch_size)):
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

                isol_event = muon_isolation_one_event(
                    muon_eta_event, muon_phi_event,
                    jTower_eta_event, jTower_phi_event, jTower_et_event,
                    threshold
                )
                #and append it to the result
                res.append(isol_event)

    return res
# %%
