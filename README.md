# Beyond Recombination: Exploring the Impact of Meiotic Frequency on Genome-wide Genetic Diversity

## Content 
**/figures** contains all the R scripts to reproduce the plots and statistical tests from the paper (use dataset from **/data**).  
**/data** contains the merged datasets obtained from all the experiments.
**/scripts** contains all the scripts used run the simulations and the individual results.  

## About the simulations 

The main script of the pipeline is **sln_model.sh** which will run SLiM (simulations) and a python script (summary statistics). It also contains all the parameters for simulation: 

- Ne: Population size
- mu: Mutation rate
- s: Selection coefficient
- h: Dominance coefficient
- n: Number of samples (for the summary statistics)
- windows: Number of windows for sumstats computation
- GR: correspond to the inverse meiotic frequency (1/GR = alpha, for GR = 1 for **alpha = 1**, 10 for **alpha = 0.01**, 50 for **alpha = 0.05** and 100 for **alpha = 0.01**)
- rho: Recombination rate (not scaled with GR)
- rep: replicate number (1 to 500)

## Usage 
### Initialization
These commands have to be run only once to setup the pipeline. It is made for SLURM based HPC.

#### Cloning the git repository
```
git clone "https://github.com/Louis-XIV-bis/meiotic-freq-simulation"
cd meiotic-freq-simulation
```

#### Create the appropriate environment using the conda export file provided
```
conda env create -f env.yml -n your_env_name
```

### Running the simulations
#### Activate the environment
```
conda activate your_env_name 
```

#### Run the same simulations than the paper 
Move to the folder of the experiement you want to reproduce (e.g. **low_s**).

Note: you may need to give the permission to some scripts (using chmod +x for instance).  

```
cd scripts/low_s
sbatch sln_model.sh
``` 

Then 500 jobs will run, one job will do the simulation for all the alpha (/ GR) values

#### Run a new simulation 
Create a folder, be careful to the name because it'll be used for the results (e.f **low_mu**).
Copy the content of the **template** folder into that new folder (you do not need the R scripts in **template/plot/**). 

Makes your changes (the parameters in sln_model.sh, changing some parameters in the .slim scripts, etc) and run the simulations.

```
cd scripts/
mkdir low_mu
cp -r template/* low_mu/
cd low_mu/
nano sln_model.sh # change mu=1e-8 to 1e-7 for instance
sbatch sln_model.sh
```
Then 500 jobs will run, one job will do the simulation for all the alpha (/ GR) values (You can also change these values in **sln_model.sh**)

### Output : results and plots

First, check that your environement is ON or change the line in **sort_files.sh** to activate it. 
When simulations are finished, run the **sort_files.sh** script to merge and removes all the intermediate results into the **/results** folder.

```
./sort_files.sh
```

For the new simulations, they will be automatically produced when running **sort_files.sh**. If you did new simulations, you'll have to add the name for the experiement (e.g. **low_mu**) at the line to load the file.

To reproduce the plots from the paper, run the Rscirpts in the **figures/** folder.

