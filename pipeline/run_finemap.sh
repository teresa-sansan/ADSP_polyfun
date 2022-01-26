#!/bin/bash

for i in {1..22} 
do
        echo "run chr$i" 
        sbatch --export=chr=$i finemap_max_snp_1.sh 
        #sbatch --export=chr=$i,max_num_snp=3 finemap_multi.sh 
        #sbatch --export=chr=$i,max_num_snp=5 finemap_multi.sh 
        #sbatch --export=chr=$i,max_num_snp=7 finemap_multi.sh 
        #sbatch --export=chr=$i,max_num_snp=10 finemap_multi.sh 

done

