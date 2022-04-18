#!/bin/bash

for max_snp in 1 3 5 7 10
do
	for i in {1..22} 
	do
		if true; then
		sbatch /gpfs/commons/home/tlin/script/polyfun/finemap_job/bellenguez/fixed_0224/bellenguez_max_snp_${max_snp}_chr${i}.sh
		fi
	done		
done

