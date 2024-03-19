#!/usr/bin/python

## Author : Louis OLLIVIER ~ louis.ollivier@etu.univ-rouen.fr 
## Université Paris-Saclay
## Lab : LISN ~ UMR9015 ~ BIOINFO team 

# This program will use the .trees files (= tree sequence, TS) given by SLiM to compute summary statistics (SFS, pi, LD, ...)

import msprime, pyslim, tskit, numpy as np, pandas as pd, sys, warnings

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

ts = tskit.load(f"trees/ts_{rep_nb}.trees") #
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

## PI ##
pi = osrts.diversity(mode='branch') # Global sumstats

# Storage of the summary statistics
pi=f"1,{pi}\n"
with open(f"pi_EGR.csv", 'a') as file:
    file.write(pi)

