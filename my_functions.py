# %%
import uproot
import numpy as np #math and science package
import scipy as sp #math and science package
import awkward as ak #root files are usuallt awkward arrays 
import matplotlib.pyplot as plt #plot stuff
import itertools

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

    #This function is designed to create a comparative plot of N sets of data with their respective errorbars included.
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
    bottom=[]
    total=[]
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

def invariant_mass(pt1, eta1, phi1, pt2, eta2, phi2):

    #This function calculates the invariant mass of the system. The first block is for muon 1 and the second for muon 2.
    #The third block adds up the contributions from muon 1 and 2. Finally, the invariant mass is computed.
    
    px1 = pt1 * np.cos(phi1)
    py1 = pt1 * np.sin(phi1)
    pz1 = pt1 * np.sinh(eta1)
    E1  = pt1 * np.cosh(eta1)

    px2 = pt2 * np.cos(phi2)
    py2 = pt2 * np.sin(phi2)
    pz2 = pt2 * np.sinh(eta2)
    E2  = pt2 * np.cosh(eta2)

    px = px1 + px2
    py = py1 + py2
    pz = pz1 + pz2
    E  = E1 + E2

    mass_squared = E**2 - (px**2 + py**2 + pz**2)
    return np.sqrt(mass_squared) if mass_squared > 0 else 0.0


def invariant_mass_all_muons(pt, eta,phi):
    #This function takes the entire set of data that contains information about the eta, phi and pt of the muons
    #and computes a list that contains the invariant mass associated to every event

    j=0

    invariant_masses=[]
    for i in range(len(eta)):
        if len(eta[i]) & len(phi[i]) & len(pt[i]) == 2:
            eta1=(eta[i])[0]
            eta2=(eta[i])[1]
            phi1=(phi[i])[0]
            phi2=(phi[i])[1]
            pt1=(pt[i])[0]
            pt2=(pt[i])[1]
            invariant_masses.append(invariant_mass(pt1, eta1, phi1, pt2, eta2, phi2))
            j=j+1
    return ak.Array(invariant_masses)

def pair_selector(pt,eta,phi,i):

 #This function aims to select the "real" pair of muons when an event has more than two muons involved

 #Let's hardcode the theoretical value of a Z mass
 ztheo=911880 #MeV

 mu=[]
 res=[]
 diff=[]
 #First, check if the lenght of the different variables is the same for some [i] event. Otherwise ignore the event, since the data
 #is not consistent 
 if (len(eta[i]) == len(phi[i]) == len(pt[i]) and len(eta[i])>2): 
    for j in range(len(eta[i])):
        aux=[(pt[i])[j],(eta[i])[j],(phi[i])[j]] #for each muon collect mu=(pt,eta,phi)
        mu.append(aux)                           #and create a vector that contains the values (mu1,mu2,mu3,...,muN)
    pairs=list(itertools.combinations(mu,2))     #create list of all possible pairs
    for mu1, mu2 in pairs:
        #for every pair take the first mu and assign pt1,eta1,phi1
        pt1, eta1, phi1 = mu1
        pt2, eta2, phi2 = mu2
        #compute the invariant mass of the pair
        a=invariant_mass(pt1, eta1, phi1, pt2, eta2, phi2)
        res.append(a)
        #and also the difference with respect to the theoretical Z mass
        b=abs(ztheo-a)
        diff.append(b)
    #Now select the best pair (meaning: pair with the smallest differencie with respect to the theoretical value)
    best_pair_position=diff.index(min(diff)) #This is the position of the best pair 
    invariant_mass_best_pair=res[best_pair_position] #The invariant mass of the best pair
    best_pair=pairs[best_pair_position] #The vector [mu1,mu2] that contains the values of pt, eta and phi of the best pair
    
    return(invariant_mass_best_pair,best_pair)

def invariant_mass_all_muons2(pt,eta,phi):

    #Improvement of invariant_mass_all_muons. This function also takes into account events with more than 2 muons. It may take much 
    #longer to run though
    #This function takes the entire set of data that contains information about the eta, phi and pt of the muons
    #and computes a list that contains the invariant mass associated to every event

    invariant_masses=[]
    
    #Scan all events
    for i in range(len(eta)):
        #If the event involves two muons
        if len(eta[i]) and len(phi[i]) and len(pt[i]) == 2:
            eta1=(eta[i])[0]
            eta2=(eta[i])[1]
            phi1=(phi[i])[0]
            phi2=(phi[i])[1]
            pt1=(pt[i])[0]
            pt2=(pt[i])[1]
            invariant_masses.append(invariant_mass(pt1, eta1, phi1, pt2, eta2, phi2))
        #If it involves more than 2 muons
        else:
            m = pair_selector(pt, eta, phi, i)
            if m is not None:  # Ensure that m is not None
                invariant_masses.append(m[0])

    return ak.Array(invariant_masses)

#-----------------------------------------------------------------------------------------------

def quality_locator(data, value):

    #This function gets the dataset that states the quality of the muons (muon_quality) and lets us 
    #select the quality using the input "value". Then it will create an array that contains a list
    #of the positions (indices) of the elements that have the desired quality

    j=0
    desired_quality_position=[]
    for index, array in enumerate(data):
        if all(x==value for x in array):
            desired_quality_position.append(index)
    return(desired_quality_position)

def quality_selector(data1, data2 ,value):

    #data1 is "muon_quality", data2 should be things like "muon_eta" or "muon_e"...
    #Using the "quality_locator" function, this function returns a list of data that involve only the 
    #muons of selected quality. For example:
    #say data2="muon_phi" and value=0, it will return a list with the values of muon_phi only for quality 0 muons
     
    selected_data=[]
    desired_quality_position=quality_locator(data1,value)
    for i in range(len(desired_quality_position)):
        pos=desired_quality_position[i]
        selected_data.append(data2[pos])
    print("Only", (len(selected_data)/len(data2))*100, r"% of the data has been selected")
    return(ak.Array(selected_data))

#-----------------------------------------------------------------------------------------------

# %%
