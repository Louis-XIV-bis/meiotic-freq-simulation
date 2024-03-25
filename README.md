# Beyond Recombination: Exploring the Impact of Meiotic Frequency on Genome-wide Genetic Diversity

## Content 
**/plot_paper** contains all the R scripts and data to reproduce the plots and statistical tests from the paper.  
**/script_simu** contains all the scripts to run the simulations and the results (more details below).  

## About the simulations 
For each tested conditions (e.g. varying s or rho), the different test **m** (meiotic frequency) values are in separate folder. **EGR** for **m = 1**, **10GR** for **m = 0.1**, **50GR** for **m = 0.05** and **100GR** for **m = 0.01**. 

sln_model.sh
This script launches SLiM simulations to generate tree sequences and then computes substats using Python.

Parameters:
Ne: Population size
mu: Mutation rate
s: Selection coefficient
h: Dominance coefficient
n: Number of samples
f: Frequency at which mutation appears under selection (in SLiM)
windows: Number of windows for sumstats computation
GR: Growth rate (10GR / 100GR / 1000GR)
rho: Recombination rate (not scaled with GR)
rep: Repetition number

## How to run 



