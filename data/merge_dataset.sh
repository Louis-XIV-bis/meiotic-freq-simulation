#!/bin/bash

# This script will gather all the data from individual experiments into a single one to do the plots for the paper

cp ../ctrl/data/pi_full_chr2_ctrl.csv . # specific for fig2

# Function to merge CSV files based on a pattern (e.g., "pi" or "t")
merge_csv_files() {
    local pattern=$1
    local exp_list=("ctrl" "h0p2" "h0p8" "high_rho" "high_s" "low_rho" "low_s" "rho_fixe")
    local output_file="${pattern}_merged.csv"
    local header_written=false

    for exp_name in "${exp_list[@]}"; do
        local results_dir="../script/${exp_name}/results"
        local csv_file="${results_dir}/${pattern}_${exp_name}.csv"

        if [ -f "$csv_file" ]; then
            if [ "$header_written" = false ]; then
                head -n 1 "$csv_file" > "$output_file"
                header_written=true
            fi
            tail -n +2 "$csv_file" >> "$output_file"
        else
            echo "File $csv_file not found!"
        fi
    done

    echo "Merging complete! Output file is: $output_file"
}

# Call the function with different patterns
# merge_csv_files "pi"
# merge_csv_files "t"
merge_csv_files "t_pichr2_full"
