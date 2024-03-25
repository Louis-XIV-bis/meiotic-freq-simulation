#!/bin/bash

# Function to replace all occurrences of "sbatch " with "./" in a file
replace_sbatch() {
    local file="$1"
    # Use sed to replace "sbatch " with "./"
    sed -i 's/sbatch /.\//g' "$file"
}

# Main function to traverse through directories and process files
traverse_directory() {
    local directory="$1"

    # Check if directory exists
    if [ -d "$directory" ]; then
        # Loop through files and directories
        for entry in "$directory"/*; do
            if [ -d "$entry" ]; then
                # Recursive call for subdirectories
                traverse_directory "$entry"
            elif [ -f "$entry" ]; then
                # Process files named "job.sh"
                if [ "$(basename "$entry")" = "jobs.sh" ]; then
                    # Replace "sbatch " with "./"
                    replace_sbatch "$entry"
                    echo "Replaced 'sbatch ' with './' in: $entry"
                fi
            fi
        done
    else
        echo "Directory '$directory' does not exist."
        exit 1
    fi
}

# Check if correct number of arguments are provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <directory>"
    exit 1
fi

# Store argument in a variable
directory="$1"

# Call the main function to traverse directories and process files
traverse_directory "$directory"
