#! /bin/bash

## Author : Louis OLLIVIER ~ louis.ollivier@etu.univ-rouen.fr 
## Université Paris-Saclay
## Lab : LISN ~ UMR9015 ~ BIOINFO team 

for rep in {1..100}
do
	#sbatch neutral_model.sh $rep	
	./neutral_model.sh $rep
done

