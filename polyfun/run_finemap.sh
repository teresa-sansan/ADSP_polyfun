#!/bin/bash 

for i in {1..22}
do
	for max_snp in  1 5 10
	do	
		echo "run_max_snp_$max_snp"
        	echo "run chr$i"

		if true; then
        	sbatch --export=chr=$i,max_num_snp=${max_snp} finemap.sh
        fi
	done
done
