#!/bin/bash
#### rerun wightman_susie
## run bellenguez susie 0817

## run jansen

for i in {1..22}
do
	for max_snp in 1 5 
	do	
		echo "run_max_snp_$max_snp"
        	echo "run chr$i"

		if false; then
		sbatch --export=chr=$i,max_num_snp=$max_snp polyfun_finemap_susie.sh 
		fi

		if true; then
		echo run jansen
        	sbatch --export=chr=$i,max_num_snp=${max_snp} finemap.sh
        	fi

	        if false; then 
		echo run kunkle
        	sbatch --export=chr=$i,max_num_snp=${max_snp} finemap_kunkle.sh 
		fi

		if false; then
        	echo run wightman
        	sbatch --export=chr=$i,max_num_snp=${max_snp} finemap_wightmansusie.sh
        fi
	done
		
done
