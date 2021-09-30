#!/bin/bash 
       
sbatch --export=max_snp=3 finemap_rescue.sh 
sbatch --export=max_snp=5 finemap_rescue.sh 
sbatch --export=max_snp=7 finemap_rescue.sh 
sbatch --export=max_snp=1 finemap_rescue.sh 



