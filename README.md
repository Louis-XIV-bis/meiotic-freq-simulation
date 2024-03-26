# Beyond Recombination: Exploring the Impact of Meiotic Frequency on Genome-wide Genetic Diversity

## Content 
**/plot_paper** contains all the R scripts and data to reproduce the plots and statistical tests from the paper.  
**/script_simu** contains all the scripts to run the simulations and the results (more details below).  

## About the simulations 
For each tested conditions (e.g. varying s or rho), the different test **m** (meiotic frequency) values are in separate folder. **EGR** for **m = 1**, **10GR** for **m = 0.1**, **50GR** for **m = 0.05** and **100GR** for **m = 0.01**. 

The main script of the pipeline is **sln_model.sh** which will run SLiM (simulations) and a python script (summary statistics). It also contains all the parameters for simulation: 

- Ne: Population size
- mu: Mutation rate
- s: Selection coefficient
- h: Dominance coefficient
- n: Number of samples (for the summary statistics)
- f: Frequency at which mutation appears under selection (in SLiM)
- windows: Number of windows for sumstats computation
- GR: meiotic frequency (1 for **m = 1**, 10 for **m = 0.01**, 50 for **m = 0.05** and 100 for **m = 0.01**)
- rho: Recombination rate (not scaled with GR)
- rep: replicate number (1 to 200)

## Usage 
### Initialization
These commands have to be run only once to setup the pipeline.

#### Cloning the git repository
```
git clone "https://github.com/Louis-XIV-bis/meiotic-freq-simulation"
cd meiotic-freq-simulation
```

#### Create the appropriate environment using the conda export file provided
```
conda env create -f scripts_simu/env.yml -n your_env_name
```

### Running the simulations
#### Activate the environment
```
conda activate your_env_name 
```
#### Run the simulations
Move to the folder of the conditions you want to test (e.g. **scripts_simu/two_chromosomes**).

Note: you may need to give the permission to some scripts (using chmod +x for instance).  

```
cd scripts_simu/two_chromosomes
./run.sh
```

This command will run consecutively the simulation for each **m** values in the terminal. 

### Output : results and plots
When finished, run the **sort_files.sh** script to merge and move all the results files a folder files within each **m** folder. 
```
./sort_files.sh
```

Additionnaly, all the results of each statistics will be copied in a specific folder within the **/plot** folder. In each folder, you'll find R script to produce basic plot of the data (not exactly the same as the paper). 
 

