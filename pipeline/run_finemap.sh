#!/bin/bash

for i in {1..22}
do
        echo "run chr$i" 
        sbatch --export=chr=$i finemap_kunkle.sh 
       # sbatch --export=chr=$i,max_num_snp=3 finemap_kunkle_multi.sh 
       # sbatch --export=chr=$i,max_num_snp=5 finemap_kunkle_multi.sh 
       # sbatch --export=chr=$i,max_num_snp=7 finemap_kunkle_multi.sh 
       # sbatch --export=chr=$i,max_num_snp=10 finemap_kunkle_multi.sh 

done

