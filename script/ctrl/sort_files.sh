#!/bin/bash

## Université Paris-Saclay
## Lab : LISN ~ UMR9015 ~ BIOINFO team

source /shared/ifbstor1/software/miniconda/bin/activate SLiMLouis

exp_name=$(basename "$PWD")

string_list=("100" "50" "10" "1")

for GR in "${string_list[@]}"; do
    
    > results/${exp_name}_chr1_${GR}.csv
    > results/${exp_name}_chr2_${GR}.csv
    > results/${exp_name}_fullchr2_${GR}.csv
          
    for rep in {1..500}; do
        cat "results/pi_chr1_${rep}_${GR}.csv" >> results/${exp_name}_chr1_${GR}.csv
        echo '' >> results/${exp_name}_chr1_${GR}.csv

        cat "results/pi_chr2_${rep}_${GR}.csv" >> results/${exp_name}_chr2_${GR}.csv
        echo '' >> results/${exp_name}_chr2_${GR}.csv
        
        cat "results/pi_fullchr2_${rep}_${GR}.csv" >> results/${exp_name}_fullchr2_${GR}.csv
        echo '' >> results/${exp_name}_fullchr2_${GR}.csv

        grep '#OUT' results/runInfo_${rep}_${GR}.txt  | cut -d' ' -f2 >> results/${exp_name}_tfix_$GR.txt;
    done
done

python pi_files.py $exp_name
python t_files.py $exp_name
Rscript ../template/plot/pi.R ${exp_name}
Rscript ../template/plot/fixation.R ${exp_name}

rm results/pi_chr1* results/pi_chr2*
rm slurm-*
rm results/runInfo_*


