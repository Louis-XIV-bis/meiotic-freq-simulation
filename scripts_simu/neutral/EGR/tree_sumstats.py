#!/usr/bin/python

## Author : Louis OLLIVIER ~ louis.ollivier@etu.univ-rouen.fr 
## Université Paris-Saclay
## Lab : LISN ~ UMR9015 ~ BIOINFO team 

# This program will use the .trees files given by SLiM to compute the SFS 
# and then do a demography parameters inference with dadi (thanks to this SFS).

import msprime, pyslim, tskit, numpy as np, pandas as pd, allel, sys

rep_nb = int(sys.argv[1])  # Replicate number
n = int(sys.argv[2])  # number of sample (diploïd) 
ne = int(sys.argv[3]) # Size of the population 
mu = float(sys.argv[4]) # Mutation rate
rho = float(sys.argv[5]) # Recombination rate

rep="EGR_"+str(rep_nb) # Replicate ID /!\ OUTPUT NAME HERE /!\
print(rep)

##### SAMPLING, ADDING NEUTRAL MUTATION AND GETTING SNP MATRIX #####

ts = pyslim.load(f"neutral_model_{rep_nb}.trees")

# Recapitate the simulation to provide a “prior history” for the initial generation of the simulation = fully coalesced
tsRecap = ts.recapitate(recombination_rate=rho, Ne=ne) # Crossing over recombination set to 0.

orig_max_roots = max(t.num_roots for t in ts.trees())

recap_max_roots = max(t.num_roots for t in tsRecap.trees())
# ~ print(f"Maximum number of roots before recapitation: {orig_max_roots}\n"
      # ~ f"After recapitation: {recap_max_roots}")

# Simplify to a subset of the population that is still alive and 
# returns a simplified tree sequence that retains only the history of the samples.

sampleInds = np.random.choice(tsRecap.individuals_alive_at(0), n, replace=False) 

# Get the node of our samples 
sampleNodes = []
for i in sampleInds:
   sampleNodes.extend(tsRecap.individual(i).nodes)
  
tsRecap_Sampled = tsRecap.simplify(sampleNodes) # This return the simplify tree (subtree) for the sampled individuals 

# Add neutral mutations (= overlay)
tsRecap_Sampled_Mutated = pyslim.SlimTreeSequence(msprime.mutate(tsRecap_Sampled, rate=(mu), keep=True)) # keep existing mutations

##### SUMMARY STATISTICS (GLOBAL) #####

## SFS ##
# Raw SFS 
rawSFS = tsRecap_Sampled_Mutated.allele_frequency_spectrum(mode="site", windows=None, polarised=True, span_normalise=False) 
# The first element (0) is not included in the SFS (see : https://tskit.dev/tutorials/analysing_tree_sequences.html#sec-tutorial-afs)
# and the last one, the fixed mutations (???)

freq_rawSFS=rawSFS[1:-1]/sum(rawSFS)
freq_rawSFS.tofile(f"rawSFS_tmp_{rep}.csv",sep=',') # In order to plot the SFS, we don't keep this element

# ~ print(rawSFS) 
# ~ print(rawSFS[1:-1]/sum(rawSFS)) 

# Normalized SFS

# Now we want to calculate the exptected SFS (under a neutral model) and then normalize the obsSFS with this value 

def expectedSFS(n): # Return the expected frequency for 1/2n to (n-1)/2n

	cumsum = 0 
	for i in range(1,2*n) :
		cumsum = cumsum + (1/i)

	exptSFS=[] # It will contain the expected frequency 
	[exptSFS.append((1/j)/cumsum) for j in range(1,2*n)] 
	
	return(exptSFS)

expSFS=np.array(expectedSFS(n),dtype=float)
# ~ print(expSFS)

normSFS = freq_rawSFS/expSFS # Calculate the normalized SFS 
# ~ print(normSFS) 

normSFS.tofile(f"normSFS_tmp_{rep}.csv",sep=',') 

## PI ##
pi = tsRecap_Sampled_Mutated.diversity()

## TAJIMA'S D ## 
D = tsRecap_Sampled_Mutated.Tajimas_D()

## WATTERSON'S THETA ## 
watterson_theta = tsRecap_Sampled_Mutated.segregating_sites()

## DAF ## 
derived_allele_freq=[i/(2*n) for i in range(1,2*n)] # Frequency 1/40 to 39/40  

DAF=0
for i in range(1,2*n-1):
	DAF+=derived_allele_freq[i]*freq_rawSFS[i]

# Storage of the summary statistics

sum_stats=np.array([pi,D,watterson_theta,DAF])
sum_stats.tofile(f"sumstats_tmp_{rep}.csv", sep=',')

##### SUMMARY STATISTICS (WINDOWS) #####

windows_sumstats = np.linspace(0, tsRecap_Sampled_Mutated.sequence_length, num=100)
# ~ print(windows_sumstats)

## PI ##
pi_windows = tsRecap_Sampled_Mutated.diversity(windows=windows_sumstats, mode='branch')
pi_windows.tofile(f"pi_tmp_{rep}.csv", sep=',')

## TAJIMA'S D ## 
D_windows = tsRecap_Sampled_Mutated.Tajimas_D(windows=windows_sumstats, mode='branch')
D_windows.tofile(f"D_tmp_{rep}.csv", sep=',')

## WATTERSON'S THETA ## 
watterson_theta_windows = tsRecap_Sampled_Mutated.segregating_sites(windows=windows_sumstats, mode='branch')
watterson_theta_windows.tofile(f"Wtheta_tmp_{rep}.csv", sep=',')

