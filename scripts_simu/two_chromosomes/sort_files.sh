# A lancer depuis le doosier contenant les differents sous dossiers (EGR, 100GR, plot etc)

string_list=("100GR" "EGR")

for GR in "${string_list[@]}"; do

	cd $GR
	
    mkdir results/trees
	mv *.trees ./results/trees

	mkdir results/runInfo
	mv runInfo* ./results/runInfo

	for i in {1..200}; do cat "D_tmp_${i}.csv" >> ./results/D_$GR.csv; echo '' >> ./results/D_$GR.csv; done
	for i in {1..200}; do cat "pi_tmp_${i}.csv" >> ./results/pi_$GR.csv; echo '' >> ./results/pi_$GR.csv; done
	for i in {1..200}; do cat "Wtheta_tmp_${i}.csv" >> ./results/Wtheta_$GR.csv; echo '' >> ./results/Wtheta_$GR.csv; done
	for i in {1..200}; do cat "normSFS_tmp_${i}.csv" >> ./results/normSFS_$GR.csv; echo '' >> ./results/normSFS_$GR.csv; done
	for i in {1..200}; do cat "rawSFS_tmp_${i}.csv" >> ./results/rawSFS_$GR.csv; echo '' >> ./results/rawSFS_$GR.csv; done
	for i in {1..200}; do cat "sumstats_tmp_${i}.csv" >> ./results/sumstats_$GR.csv; echo '' >> ./results/sumstats_$GR.csv; done
	for i in {1..200}; do cat "rawSFS_local_tmp_${i}.csv" >> ./results/rawSFS_local_$GR.csv; echo '' >> ./results/rawSFS_local_$GR.csv; done
	for i in {1..200}; do cat "b_one_tmp_${i}.csv" >> ./results/b_one_$GR.csv; echo '' >> ./results/b_one_$GR.csv; done
	for i in {1..200}; do cat "sackin_tmp_${i}.csv" >> ./results/sackin_$GR.csv; echo '' >> ./results/sackin_$GR.csv; done
	
	rm rawSFS_tmp_* normSFS_tmp_* sumstats_tmp_* D_* pi_* Wtheta_* rawSFS_local_* sackin_* b_one_*
	
	cd results
	cp * ../../plot # pas -r donc ne prend pas les dossiers (runInfo etc) 

	cd ../.. # retour au dossier de base 
done

cd ./plot/

mv D* sumstats/windows/D
mv Wtheta* sumstats/windows/Wtheta
mv pi* sumstats/windows/pi
mv b_one* sumstats/windows/tree_balance
mv sackin* sumstats/windows/tree_balance
mv sumstats_* sumstats/nowindow
mv rawSFS_local* SFS/raw_local
mv rawSFS* SFS/raw
mv normSFS* SFS/normalized

