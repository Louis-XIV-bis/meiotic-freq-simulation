string_list=("100GR" "EGR")

for GR in "${string_list[@]}"; do
	for i in {1..200}; do
		grep '#OUT' ../$GR/results/runInfo/runInfo_$i.txt  | cut -d' ' -f2 >> fix_$GR.txt;
	done
done
