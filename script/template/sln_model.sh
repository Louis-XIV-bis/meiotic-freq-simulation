#!/bin/bash

## Université Paris-Saclay
## Lab : LISN ~ UMR9015 ~ BIOINFO team

#SBATCH --job-name=template
#SBATCH --array=1-500
#SBATCH --mem-per-cpu=32GB
#SBATCH --account=yeast_neutral_model

source /shared/ifbstor1/software/miniconda/bin/activate SLiMLouis

## This script will launch the SLiM simulations that will generate tree sequence and then 
## use Python to calculate some sumstats for the given model.

# Parameters for the simulations

Ne=1000 # Population size 
mu=1e-8 # Mutation rate 
s=0.05 # Selection coefficient
h=0.5 # Dominance coefficient 
n=100 # Number of samples for the python part
windows=50 # Number of windows for sumstats computation on the first chromosome
rho=5e-8 # Recombination rate (not scaled with GR!)

rep=${SLURM_ARRAY_TASK_ID} # Replicate number

# Test different value for GR (= 1/alpha) 
for GR in 100 50 10 1; do 

    # Run the SLiM simulation : generate a .trees file (= a tree sequence)
    tmpSLiM="results/runInfo_${rep}_${GR}.txt" # Will contains informations about the run (seed, time, etc)

    fixed=0

    while ((fixed==0))
    do 
        slim -d Ne=$Ne -d rho=$rho -d rep=$rep -d s=$s -d h=$h -d mu=$mu -d GR=$GR sln_sweep.slim > $tmpSLiM
        fixed=$(grep m2 $tmpSLiM|wc -l) 
        echo $fixed
    done

    # Run the Python script : compute SFS, do the demography parameters inference and save them in a file
    python3 tree_sumstats.py $rep $n $Ne $mu $rho $windows $GR

    rm results/ts_${rep}_${GR}.trees
done
