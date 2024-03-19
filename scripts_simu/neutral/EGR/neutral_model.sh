#!/bin/bash

#SBATCH --job-name=neutral_EGR
#SBATCH --mem=128GB
#SBATCH -p long

source /shared/ifbstor1/software/miniconda/bin/activate SLiMLouis

## Author : Louis OLLIVIER ~ louis.ollivier@etu.univ-rouen.fr
## Université Paris-Saclay
## Lab : LISN ~ UMR9015 ~ BIOINFO team

## This script will launch the SLiM simulations that will generate tree sequence and then 
## use Python to calculate some substats for the given model.

# Parameters for the simulations

Ne=1000 # Population size 
mu=1e-8 # Mutation rate 
rho=5e-8 # Recombination rate 
# s=0.01 # Selection coefficient
n=20 # Number of samples

rep=$1 

# Run the SLiM simulation : generate a .trees file (= a tree sequence)

tmpSLiM="runInfo_${rep}.txt" # Will contains informations about the run (seed, time, etc)

slim -d Ne=$Ne -d rho=$rho -d rep=$rep neutral_EGR.slim > $tmpSLiM

# seed=$(sed -n "2p" $tmpSLiM) # The seed will be available in case of problem
# echo $seed > "seed_${rep}.txt"

# Run the Python script : compute SFS, do the demography parameters inference and save them in a file

python3 tree_sumstats.py $rep $n $Ne $mu $rho

# Delete all the temporary files 
# rm neutral_model_${rep}.trees 
rm runInfo_${rep}.txt # Contains informations about memory usage, time etc if needed
