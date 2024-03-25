cd 100GR/
./jobs.sh

cd ../EGR/
./jobs.sh
./ts.sh "#SBATCH --mem=128GB"
./ts.sh "source /shared/ifbstor1/software/miniconda/bin/activate SLiMLouis"
./ts.sh "## Author : Louis OLLIVIER ~ louis.ollivier@etu.univ-rouen.fr"



