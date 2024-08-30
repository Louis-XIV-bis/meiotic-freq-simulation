#!/usr/bin/python

## Université Paris-Saclay
## Lab : LISN ~ UMR9015 ~ BIOINFO team 

# This program will use the .trees files (= tree sequence, TS) given by SLiM to compute summary statistics (pi)

import msprime, pyslim, tskit, numpy as np, pandas as pd, sys, warnings

rep = int(sys.argv[1]) # Replicate number
n = int(sys.argv[2]) # Number of sample for the simplication and sumstats (diploïd) 
ne = int(sys.argv[3]) # Size of the population 
mu = float(sys.argv[4]) # Mutation rate
rho = float(sys.argv[5]) # Recombination rate
nb_windows = int(sys.argv[6]) 
GR = int(sys.argv[7])

##### SAMPLING, ADDING NEUTRAL MUTATION AND GETTING SNP MATRIX #####
##### RECAPITATE, SAMPLING AND ADDING NEUTRAL MUTATION) ###################################
# https://tskit.dev/pyslim/docs/latest/tutorial.html

ts = tskit.load(f"results/ts_{rep}_{GR}.trees")
ts = pyslim.update(ts)

warnings.simplefilter('ignore', msprime.TimeUnitsMismatchWarning)

# Recapitate the simulation to provide a “prior history” for the initial generation of the simulation = fully coalesced #
# RTS : recapitated 

rts = pyslim.recapitate(ts, recombination_rate=rho, ancestral_Ne=ne)
      
# Sampling indivuals from the TS to compute sumstats later # 
# SRTS = simplified recapitated TS 

rng = np.random.default_rng()
alive_inds = pyslim.individuals_alive_at(rts, 0) # sample individuals alive = coming from the last generation
keep_indivs = rng.choice(alive_inds, n, replace=False)
keep_nodes = []
for i in keep_indivs:
  keep_nodes.extend(rts.individual(i).nodes)
  
srts = rts.simplify(keep_nodes, keep_input_roots=True)

# Adding neutral mutations to the samples (overlay) # 
# OSRTS = overlay simplified recapitated TS 

osrts = msprime.mutate(srts, rate=(mu), keep=True) # keep existing mutations

##### SUMMARY STATISTICS (WITH WINDOWS) ####################################

# We want a user defined number of windows for the first chromosome but only 1 value
# for the second since it's entirely neutral 

windows_chr = np.linspace(0, osrts.sequence_length/2, num=nb_windows)
windows_chr = [int(x) for x in windows_chr] + [int(osrts.sequence_length)]

## PI ##
pi = osrts.diversity(windows=windows_chr, mode='site')
pi_chr1 = pi[:-1]
pi_chr2 = pi[-1]

pi_chr1.tofile(f"results/pi_chr1_{rep}_{GR}.csv", sep=',')
pi_chr2.tofile(f"results/pi_chr2_{rep}_{GR}.csv", sep=',')
