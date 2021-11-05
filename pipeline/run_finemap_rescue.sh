#!/bin/bash 

#sbatch finemap_rescue_snp1_all_2.sh  
sbatch --export=max_snp=3 finemap_rescue.sh 
sbatch --export=max_snp=5 finemap_rescue.sh 
sbatch --export=max_snp=7 finemap_rescue.sh 
sbatch --export=max_snp=10 finemap_rescue.sh 



