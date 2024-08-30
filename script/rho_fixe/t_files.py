#!/usr/bin/python

## Université Paris-Saclay
## Lab : LISN ~ UMR9015 ~ BIOINFO team 

import numpy as np, re, pandas as pd, sys

# We want to gather all the data for the experiment in one file in a table with the simulation parameter and the pi 
# We'll do the mean + sd for each pi windows among the replicate for a given GR (1/alpha)

# Retrieve the value for the parameter 
variables = {}
with open('sln_model.sh', 'r') as file:
    for line in file:
        # Match lines that define variables in the script 
        match = re.match(r'^(\w+)\s*=\s*([^\s#]+)', line)
        if match:
            var_name = match.group(1)
            var_value = match.group(2)
            variables[var_name] = var_value

mu = float(variables.get('mu'))
s = float(variables.get('s'))
h = float(variables.get('h'))
rho = float(variables.get('rho'))
GR_values = [100, 50, 10, 1]

exp_name = str(sys.argv[1]) # experiment name to open the data / save the final file 

# Initiailize the data_table 
# Prepare a list to store all rows of the table
table_rows = []

# Loop through the list of GR values
for GR in GR_values:

    filename = f"results/{exp_name}_tfix_{GR}.txt"  

    try:
        df = pd.read_csv(filename)
        df = df - 2000 # neutral generations

        # Calculate mean and standard deviation (assuming only one column)
        mean_value = df.mean().iloc[0]  # Mean of the first (and only) column
        sd_value = df.std().iloc[0]    # Standard deviation of the first (and only) column
        
        # Create a row with all the values, including mean and sd with window set to 'chr2'
        row = {
            'GR': GR,
            's': s,
            'h': h,
            'rho': rho,
            'rho_scaled': rho, # WARNING : should be rho * GR in other cases
            'mean': mean_value,
            'sd': sd_value
        }
        table_rows.append(row)
    
    except FileNotFoundError:
        print(f"File {filename} not found. Skipping GR = {GR}.")

# Create a DataFrame from the list of rows
df = pd.DataFrame(table_rows)
df.to_csv(f'results/t_{exp_name}.csv', index = False)
