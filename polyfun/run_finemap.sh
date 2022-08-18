#!/bin/bash
#### rerun wightman_susie
## run bellenguez susie 0817

## run chr19-22

for i in {3..22}
do
	for max_snp in 10
	do	
		echo "run_max_snp_$max_snp"
        	echo "run chr$i"

		if false; then
		sbatch --export=chr=$i,max_num_snp=$max_snp polyfun_finemap_susie.sh 
		fi

		if true; then
		echo run bellenguez
        	sbatch --export=chr=$i,max_num_snp=${max_snp} finemap.sh
        	fi

	       if false; then 
		echo run kunkle
        	sbatch --export=chr=$i,max_num_snp=${max_snp} finemap_kunkle.sh 
		fi

		if false; then
        	echo run wightman
        	sbatch --export=chr=$i,max_num_snp=${max_snp} finemap.sh
        fi
	done
		
done
