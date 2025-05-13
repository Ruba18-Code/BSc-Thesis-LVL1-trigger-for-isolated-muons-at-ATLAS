# %%
import uproot
import numpy as np #math and science package
import scipy as sp #math and science package
import awkward as ak #root files are usuallt awkward arrays 
import matplotlib.pyplot as plt #plot stuff

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
#-----------------------------------------------------------------------------------------------
def number_mu_per_event(data):

    element_len=np.zeros(len(data), dtype=int)
    for i in range(len(data)):
        aux=len(data[i])
        element_len[i]=aux    
    event=np.zeros(np.max(element_len)+1, dtype=int)

    for i in range(np.min(element_len),np.max(element_len)):
        for j in range(len(element_len)):
            if i ==element_len[j]:
               event[i]=event[i]+1
    return(event)


def plot_number_mu_per_event(data):
    
    #The following function takes the data, then counts the lenght of every element, and stores it in the vector element_len
    #then, the second loop counts how many times does each lenght appear on element_len and creates the event vector
    #the element n of this vectors contains an integer that represents how many times does the number n appear in element_len
    #finally, we plot it 

    event=number_mu_per_event(data)

    #-------PLOT-----------------------------
    #            
    x= np.arange(len(event)) #create the x axis

    plt.bar(x, event,color=plt.cm.viridis(0.9), edgecolor= 'black', alpha=0.7) #plot
    for i in range(len(x)):
        plt.text(x[i], event[i] + 0.2, str(event[i]), ha='center', va='bottom', fontsize=8)
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.xlabel("Number of Muons recorded")
    plt.ylabel("Counts")
    plt.title("Number of recorded muons per event")
    plt.show()

def plot_number_mu_per_event2(data1,data2,label1,label2):
    event=number_mu_per_event(data1)
    x= np.arange(len(event)) #create the x axis

    plt.bar(x, event,color=plt.cm.viridis(0.9), edgecolor= 'black', alpha=0.7) #plot
    for i in range(len(x)):
        plt.text(x[i], event[i] + 0.2, str(event[i]), ha='center', va='bottom', fontsize=8)
    
    event=number_mu_per_event(data2)
    
    x= np.arange(len(event)) #create the x axis

    plt.bar(x, event,color=plt.cm.viridis(0.1), edgecolor= 'black', alpha=0.7) #plot
    for i in range(len(x)):
        plt.text(x[i], event[i] + 0.2, str(event[i]), ha='center', va='bottom', fontsize=8)

    plt.grid(True, linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.xlabel("Number of Muons recorded")
    plt.ylabel("Counts")
    plt.legend(labels=[label1,label2])
    plt.title("Number of recorded muons per event")
    plt.show()

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

def high_quality_histogram2(data1,data1_quality,data2,data2_quality,nbins,x_range,x_label,y_label,title_label,label1,label2):
    value=0
    selected_muons1=quality_selector(data1_quality,data1,value)
    selected_muons2=quality_selector(data2_quality,data2,value)

    histogram2errors(selected_muons1,selected_muons2,nbins,x_range,x_label,y_label,title_label,label1,label2)
#-----------------------------------------------------------------------------------------------

def plot(x,bins,label=None,selection=np.array([])):
    if len(selection) == 0:
        selection=np.full(len(x),True)
        
    bins_overflow = np.concatenate([[-np.inf],bins[1:-1],[np.inf]])
    hists = [np.histogram(x[(process==i)&selection],bins=bins_overflow,weights=weights[(process==i)&selection])[0]
             for i in range(10)]
    total = sum(hists)
    
    bin_index = np.digitize(x,bins_overflow)-1
    uncerts = [np.array([np.sqrt(sum(weights[(bin_index==j)&(process==i)*selection]**2)) for j in range(len(bins)-1)])
               for i in range(10)]
    total_uncerts = np.sqrt(sum([u**2 for u in uncerts]))
    
    data = np.histogram(x[isData&selection],bins=bins_overflow)[0]
    data_uncerts = np.sqrt(data)
    ratio = data/total
    ratio_uncerts = np.sqrt((data_uncerts/total)**2+(data/total**2*total_uncerts)**2)
    
    fig, ax, rax = ampl.ratio_axes()

    for i in range(10):
        ax.bar(bins[:-1],
               hists[i],
               color=colors[i],
               bottom=np.zeros(len(bins)-1) if i==0 else sum(hists[:i]),
               align='edge',
               width=bins[1:]-bins[:-1],
               edgecolor='black',
               label=labels[i])

    ax.bar(bins[:-1],
           2*total_uncerts,
           bottom=total-total_uncerts,
           color=(0,0,0,0),
           hatch='////',
           align='edge',
           width=bins[1:]-bins[:-1],label='total MC stat. uncert.')
    
    ax.errorbar(bins[:-1]+0.5*(bins[1:]-bins[:-1]),data,yerr=data_uncerts,fmt='.',color='black')

    rax.errorbar(bins[:-1]+0.5*(bins[1:]-bins[:-1]),ratio,yerr=ratio_uncerts,fmt='.',color='black')
    rax.axhline(1,linestyle='--',color='black')

    ax.set_xlim(bins[0],bins[-1])
    ax.set_ylim(0)

    if label != None:
        rax.set_xlabel(label)
    rax.set_ylabel('Data/MC')
    ax.set_ylabel('Events')
    ax.legend(frameon=False,bbox_to_anchor=(1.05,1),loc='upper left',fontsize='small')
    plt.tight_layout()


def ak_element_lenght_counter(data):
    return(ak.num(data))

def coolplot(data,bins,colors,labels,x_label,y_label,title):

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

    #Compute the histogram, error and bottom (used for the errorbars) for every dataset

    for i in range(len(data)):
        hists.append(np.histogram(ak.flatten(data[i], axis=None),bins=bins_overflow)[0])
        error.append(np.sqrt(hists[i]))
        bottom.append(hists[i]-error[i])

    #Plot according to some preferences


    for i in range(len(data)):

     error[i]=error[i]/np.sum(hists[i])
     bottom[i]=bottom[i]/np.sum(hists[i])
     hists[i]=hists[i]/np.sum(hists[i])

     plt.step(
        bins[:-1],
        hists[i],
        where='post',  # or 'mid', or 'pre' depending on visual preference
        color=colors[i],
        label=labels[i],
        alpha=0.5)
        
        #This part plots the errorbars, bottom ensures that the bars start at y=hist-error and, since bars are 2*error long,
        #half of the bar will be below the central value and half will be above it

     plt.bar(bins[:-1],
            2*error[i],
            bottom=bottom[i],
            hatch='////',
            align='edge',
            width=bins[1:]-bins[:-1],
            color=(0,0,0,0),
            linestyle='--')
    
    #cool parameters + plot the title and labels

    plt.tight_layout()
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(title)

    return(hists_norm)


# %%
