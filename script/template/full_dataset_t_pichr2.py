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
windows = int(variables.get('windows'))
rho = float(variables.get('rho'))
GR_values = [100, 50, 10, 1]
nb_rep = 500 

exp_name = str(sys.argv[1]) # experiment name to open the data / save the final file 

# Initiailize the data_table 
# Prepare a list to store all rows of the table
table_rows = []

# Loop through the list of GR values
for GR in GR_values:

    # filename_chr1 = f"results/{exp_name}_chr1_{GR}.csv" 
    filename_chr2_pi = f"results/{exp_name}_chr2_{GR}.csv"
    filename_t = f"results/{exp_name}_tfix_{GR}.txt" 
 
    try:
        df_pi = pd.read_csv(filename_chr2_pi)
        df_t = pd.read_csv(filename_t)

        for rep in range(1, nb_rep):
            pi = df_pi.iloc[rep - 1, 0] # Store the pi value for the replicate 
            t = df_t.iloc[rep - 1, 0] # and associated t

            # Create a row with all the values, including mean and sd with window set to 'chr2'
            row = {
                'GR': GR,
                's': s,
                'h': h,
                'rho': rho,
                'rho_scaled': rho * GR,
                'window': 'chr2',
                'rep': rep,
                'pi': pi,
                't': t
            }
            table_rows.append(row)
        
    except FileNotFoundError:
        print(f"File {filename_chr2_pi} or {filename_t} not found. Skipping GR = {GR}.")

# # Create a DataFrame from the list of rows
df = pd.DataFrame(table_rows)
df.to_csv(f'results/t_pichr2_full_{exp_name}.csv', index = False)