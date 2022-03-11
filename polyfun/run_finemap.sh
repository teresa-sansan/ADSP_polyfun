#!/bin/bash

for i in {1..22} 
do
        echo "run chr$i"
	if true; then
	echo run bellenguez
        sbatch --export=chr=$i finemap_bellenguez.sh
        sbatch --export=chr=$i,max_num_snp=3 finemap_bellenguez_multi.sh
        sbatch --export=chr=$i,max_num_snp=5 finemap_bellenguez_multi.sh
        sbatch --export=chr=$i,max_num_snp=7 finemap_bellenguez_multi.sh
        sbatch --export=chr=$i,max_num_snp=10 finemap_bellenguez_multi.sh
        fi


        if false; then 
	echo run kunkle
        sbatch --export=chr=$i finemap_kunkle.sh 
        sbatch --export=chr=$i,max_num_snp=3 finemap_kunkle_multi.sh 
        sbatch --export=chr=$i,max_num_snp=5 finemap_kunkle_multi.sh 
        sbatch --export=chr=$i,max_num_snp=7 finemap_kunkle_multi.sh 
        sbatch --export=chr=$i,max_num_snp=10 finemap_kunkle_multi.sh 
	fi

	if true; then
        echo run wightman
        sbatch --export=chr=$i finemap_max_snp_1.sh
        sbatch --export=chr=$i,max_num_snp=3 finemap_multi.sh
        sbatch --export=chr=$i,max_num_snp=5 finemap_multi.sh
        sbatch --export=chr=$i,max_num_snp=7 finemap_multi.sh
        sbatch --export=chr=$i,max_num_snp=10 finemap_multi.sh
        fi
		
done

