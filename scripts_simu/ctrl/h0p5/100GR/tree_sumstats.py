#!/usr/bin/python

## Université Paris-Saclay
## Lab : LISN ~ UMR9015 ~ BIOINFO team 

# This program will use the .trees files (= tree sequence, TS) given by SLiM to compute summary statistics (SFS, pi, LD, ...)

import msprime, pyslim, tskit, numpy as np, pandas as pd, allel, sys, warnings

rep_nb = int(sys.argv[1])  # Replicate number
n = int(sys.argv[2])  # number of sample (diploïd) 
ne = int(sys.argv[3]) # Size of the population 
mu = float(sys.argv[4]) # Mutation rate
rho = float(sys.argv[5]) # Recombination rate
nb_windows = int(sys.argv[6])
GR = int(sys.argv[7])

##### SAMPLING, ADDING NEUTRAL MUTATION AND GETTING SNP MATRIX #####
##### RECAPITATE, SAMPLING AND ADDING NEUTRAL MUTATION) ###################################
# https://tskit.dev/pyslim/docs/latest/tutorial.html

ts = tskit.load(f"ts_{rep_nb}.trees") #
ts = pyslim.update(ts)

warnings.simplefilter('ignore', msprime.TimeUnitsMismatchWarning)

# Recapitate the simulation to provide a “prior history” for the initial generation of the simulation = fully coalesced #
# RTS : recapitated 
print("before recap")
rts = pyslim.recapitate(ts, recombination_rate=rho*GR, ancestral_Ne=ne)
print("test1")

# ~ orig_max_roots = max(t.num_roots for t in ts.trees())
# ~ recap_max_roots = max(t.num_roots for t in tsRecap.trees())
# ~ print(f"Maximum number of roots before recapitation: {orig_max_roots}\n"
      # ~ f"After recapitation: {recap_max_roots}")
      
# Sampling indivuals from the TS to compute sumstats later # 
# SRTS = simplified recapitated TS 

rng = np.random.default_rng()
alive_inds = pyslim.individuals_alive_at(rts, 0) # sample individuals alive = coming from the last generation
keep_indivs = rng.choice(alive_inds, n, replace=False)
keep_nodes = []
for i in keep_indivs:
  keep_nodes.extend(rts.individual(i).nodes)
  
srts = rts.simplify(keep_nodes, keep_input_roots=True)

print("test2")

# ~ print(f"Before, there were {rts.num_samples} sample nodes (and {rts.num_individuals} individuals)\n"
      # ~ f"in the tree sequence, and now there are {srts.num_samples} sample nodes\n"
      # ~ f"(and {srts.num_individuals} individuals).")

# Adding neutral mutations to the samples (overlay) # 
# OSRTS = overlay simplified recapitated TS 

osrts = msprime.mutate(srts, rate=(mu), keep=True) # keep existing mutations
print("test3")

# ~ print(f"Before overlay : {osrts.num_mutations} mutations, after ovelay : {osrts.num_mutations} mutations")

##### SUMMARY STATISTICS (NO WINDOWS) ####################################

## SFS ##
# Raw SFS #
# https://tskit.dev/tutorials/analysing_tree_sequences.html#sec-tutorial-afs)

rawSFS = osrts.allele_frequency_spectrum(mode="site", windows=None, polarised=True, span_normalise=False) 

freq_rawSFS=rawSFS[1:-1]/sum(rawSFS) # We want frequencies instead of value
freq_rawSFS.tofile(f"rawSFS_tmp_{rep_nb}.csv",sep=',') 

# Normalized SFS #
# To normalize, we need to calculate the exptected SFS (under a neutral model) and then normalize the rawSFS with this value 

def expectedSFS(n): # Return the expected frequency for 1/2n to (n-1)/2n

	cumsum = 0 
	for i in range(1,2*n) :
		cumsum = cumsum + (1/i)

	exptSFS=[] # It will contain the expected frequency 
	[exptSFS.append((1/j)/cumsum) for j in range(1,2*n)] 
	
	return(exptSFS)

expSFS=np.array(expectedSFS(n),dtype=float)

normSFS = freq_rawSFS/expSFS # Calculate the normalized SFS 

normSFS.tofile(f"normSFS_tmp_{rep_nb}.csv",sep=',') 

## PI ##
pi = osrts.diversity() # Global sumstats

## TAJIMA'S D ## 
D = osrts.Tajimas_D() # Global sumstats

## WATTERSON'S THETA ## 
watterson_theta = osrts.segregating_sites() # Global sumstats

## DAF ## 
derived_allele_freq=[i/(2*n) for i in range(1,2*n)] # Frequency 1/40 to 39/40  

DAF=0
for i in range(1,2*n-1):
	DAF+=derived_allele_freq[i]*freq_rawSFS[i]

# Storage of the summary statistics
sum_stats=np.array([pi,D,watterson_theta,DAF])
sum_stats.tofile(f"sumstats_tmp_{rep_nb}.csv", sep=',')

##### SUMMARY STATISTICS (WITH WINDOWS) ####################################

windows_sumstats = np.linspace(0, osrts.sequence_length, num=nb_windows)
windows_sumstats = [int(x) for x in windows_sumstats] 

## PI ##
pi_windows = osrts.diversity(windows=windows_sumstats, mode='branch')
pi_windows.tofile(f"pi_tmp_{rep_nb}.csv", sep=',')

## TAJIMA'S D ## 
D_windows = osrts.Tajimas_D(windows=windows_sumstats, mode='branch')
D_windows.tofile(f"D_tmp_{rep_nb}.csv", sep=',')

## WATTERSON'S THETA ## 
watterson_theta_windows = osrts.segregating_sites(windows=windows_sumstats, mode='branch')
watterson_theta_windows.tofile(f"Wtheta_tmp_{rep_nb}.csv", sep=',')


## SFS  local ##
rawSFS_local = osrts.allele_frequency_spectrum(mode="branch", windows=windows_sumstats, polarised=True, span_normalise=False) 
rawSFS_local.tofile(f"rawSFS_local_tmp_{rep_nb}.csv",sep=',') 

## TREE BALANCE INDICES ## 
list_sackin = list_b1 = np.array([])

for tree in osrts.trees():
	list_sackin = np.append(list_sackin,tree.sackin_index())
	list_b1 = np.append(list_b1,tree.b1_index())

list_sackin.tofile(f"sackin_tmp_{rep_nb}.csv", sep=',')
list_b1.tofile(f"b_one_tmp_{rep_nb}.csv", sep=',')

## LD & SNP MATRIX ## 
# https://scikit-allel.readthedocs.io/en/stable/stats/ld.html
# Formating SNP matrix for LD computation 

SNP_matrix = osrts.genotype_matrix().T # individuals in rows and sites in columns
SNP_matrix = pd.DataFrame(data=SNP_matrix.astype(int))

SNP_pos = np.round(osrts.tables.asdict()["sites"]["position"]).astype(int)

SNP_matrix.to_csv(f"SNP_matrix_{rep_nb}.csv", sep=' ', header=SNP_pos, index=False)

SNP_alt = np.empty((0,len(SNP_pos)), int) # init
split_SNP_matrix = np.array_split(SNP_matrix, 2*n) # Each element will contain a row of SNP_matrix = 1 copy of the genome

# Sum of rows 2 by 2 : sum of the 2 copy of a genome in order to get biallelic variants
for i in range(0,2*n,2): # diploïd => 2*n copy ; row n°0 = copy n°1 of ind 1, row n°1 = copy n°2 of ind 1, row n°2 = copy n°1 of ind 2, ...
    tmp_alt = np.sum((split_SNP_matrix[i],split_SNP_matrix[i+1]),axis=0)
    SNP_alt = np.append(SNP_alt, tmp_alt, axis=0)

LD_r = allel.windowed_r_squared(SNP_pos, SNP_alt.T, start = SNP_pos[0], stop = SNP_pos[-1], size = 50000)

np.savetxt(f"LD_r2_values_{rep_nb}.csv", LD_r[0])
np.savetxt(f"LD_r2_windows_{rep_nb}.csv", LD_r[1])
