#!/bin/bash

for i in {1..22}
do
	for max_snp in 1 3 5 7 10
	do	
		echo "run_max_snp_$max_snp"
        	echo "run chr$i"

		if false; then
		sbatch --export=chr=$i,max_num_snp=$max_snp polyfun_finemap_susie.sh 
		fi

		if false; then
		echo run bellenguez
        	sbatch --export=chr=$i finemap_bellenguez.sh
        	sbatch --export=chr=$i,max_num_snp=3 finemap_bellenguez_multi.sh
        	sbatch --export=chr=$i,max_num_snp=5 finemap_bellenguez_multi.sh
        	sbatch --export=chr=$i,max_num_snp=7 finemap_bellenguez_multi.sh
        	sbatch --export=chr=$i,max_num_snp=10 finemap_bellenguez_multi.sh
        	fi


	        if true; then 
		echo run kunkle
        	sbatch --export=chr=$i,max_num_snp=${max_snp} finemap_kunkle.sh 
		fi

		if false; then
        	echo run wightman
        	sbatch --export=chr=$i finemap_max_snp_1.sh
        	sbatch --export=chr=$i,max_num_snp=3 finemap_multi.sh
        	sbatch --export=chr=$i,max_num_snp=5 finemap_multi.sh
        	sbatch --export=chr=$i,max_num_snp=7 finemap_multi.sh
        	sbatch --export=chr=$i,max_num_snp=10 finemap_multi.sh
        	fi
	done
		
done

