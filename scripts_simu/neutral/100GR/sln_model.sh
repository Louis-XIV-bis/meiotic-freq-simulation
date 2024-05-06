#!/bin/bash

#SBATCH --mem=128GB

source /shared/ifbstor1/software/miniconda/bin/activate SLiMLouis

## Université Paris-Saclay
## Lab : LISN ~ UMR9015 ~ BIOINFO team

## This script will launch the SLiM simulations that will generate tree sequence and then 
## use Python to calculate some substats for the given model.

# Parameters for the simulations

Ne=1000 # Population size 
mu=1e-8 # Mutation rate 
# s=0.1 # Selection coefficient
h=0.4 # Dominance coefficient 
n=20 # Number of samples
f=0.1 # Frequency at which appears the mutation under selection (slim)
windows=500 # Number of windows for sumstats computation
GR=100 # 10GR / 100GR / 1000GR 
rho=5e-8 # Recombination rate (not scaled with GR!)

rep=$1 

# Run the SLiM simulation : generate a .trees file (= a tree sequence)

tmpSLiM="runInfo_${rep}.txt" # Will contains informations about the run (seed, time, etc)

#fixed=0

#while ((fixed==0))
#do 
slim -d Ne=$Ne -d rho=$rho -d rep=$rep -d h=$h -d mu=$mu -d f=$f -d GR=$GR sln_sweep.slim > $tmpSLiM
#fixed=$(grep m2 $tmpSLiM|wc -l) 
#echo $fixed
#done

# seed=$(sed -n "2p" $tmpSLiM) # The seed will be available in case of problem
# echo $seed > "seed_${rep}.txt"

# Run the Python script : compute SFS, do the demography parameters inference and save them in a file

python3 tree_sumstats.py $rep $n $Ne $mu $rho $windows $GR

