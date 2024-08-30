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

# Create the windows for each pi (the same than the other python script)
windows_chr1 = np.linspace(0, 1000000, num = int(windows - 1))
windows_chr1 = [int(x) for x in windows_chr1]

exp_name = str(sys.argv[1]) # experiment name to open the data / save the final file 

# Initiailize the data_table 
# Prepare a list to store all rows of the table
table_rows = []

# Loop through the list of GR values
for GR in GR_values:

    filename_chr1 = f"results/{exp_name}_chr1_{GR}.csv" 

    try:
        df = pd.read_csv(filename_chr1)
        
        # Calculate mean and standard deviation by column
        mean_values = df.mean()
        sd_values = df.std()
        
        # Loop through each window in windows_chr1 with index
        for idx, window_value in enumerate(windows_chr1):
            # Create a row with all the values, including mean and sd at the index
            row = {
                'GR': GR,
                's': s,
                'h': h,
                'rho': rho,
                'rho_scaled': rho * GR,
                'window': window_value,
                'mean': mean_values.iloc[idx] if idx < len(mean_values) else None,  # Handle out-of-range index
                'sd': sd_values.iloc[idx] if idx < len(sd_values) else None  # Handle out-of-range index
            }
            table_rows.append(row)

    except FileNotFoundError:
        print(f"File {filename_chr2} not found. Skipping GR = {GR}.")

    filename_chr2 = f"results/{exp_name}_chr2_{GR}.csv" 
    
    try:
        df = pd.read_csv(filename_chr2)
        
        # Calculate mean and standard deviation (assuming only one column)
        mean_value = df.mean().iloc[0]  # Mean of the first (and only) column
        sd_value = df.std().iloc[0]    # Standard deviation of the first (and only) column
        
        # Create a row with all the values, including mean and sd with window set to 'chr2'
        row = {
            'GR': GR,
            's': s,
            'h': h,
            'rho': rho,
            'rho_scaled': rho * GR,
            'window': 'chr2',
            'mean': mean_value,
            'sd': sd_value
        }
        table_rows.append(row)
    
    except FileNotFoundError:
        print(f"File {filename_chr2} not found. Skipping GR = {GR}.")

# Create a DataFrame from the list of rows
df = pd.DataFrame(table_rows)
df.to_csv(f'results/pi_{exp_name}.csv', index = False)

# Full chr2 pi (specific to control, for fig2) 

# Create the windows for each pi (the same than the other python script)
windows_chr2 = np.linspace(1000000, 2000000, num = int(windows - 1))
windows_chr2 = [int(x) for x in windows_chr2]


# Initiailize the data_table 
# Prepare a list to store all rows of the table
table_rows = []

# Loop through the list of GR values
for GR in GR_values:

    filename_fullchr2 = f"results/{exp_name}_fullchr2_{GR}.csv" 

    try:
        df = pd.read_csv(filename_fullchr2)
        
        # Calculate mean and standard deviation by column
        mean_values = df.mean()
        sd_values = df.std()
        
        # Loop through each window in windows_chr1 with index
        for idx, window_value in enumerate(windows_chr2):
            # Create a row with all the values, including mean and sd at the index
            row = {
                'GR': GR,
                's': s,
                'h': h,
                'rho': rho,
                'rho_scaled': rho * GR,
                'window': window_value,
                'mean': mean_values.iloc[idx] if idx < len(mean_values) else None,  # Handle out-of-range index
                'sd': sd_values.iloc[idx] if idx < len(sd_values) else None  # Handle out-of-range index
            }
            table_rows.append(row)

    except FileNotFoundError:
        print(f"File {filename_fullchr2} not found. Skipping GR = {GR}.")
        
# Create a DataFrame from the list of rows
df = pd.DataFrame(table_rows)
df.to_csv(f'results/pi_full_chr2_{exp_name}.csv', index = False)
        
