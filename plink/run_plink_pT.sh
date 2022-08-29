#!/bin/bash

for i in {1..21}
do
  echo "start calculating PRS for chr $i"
  sbatch --export=chr=$i plink_pT.sh   
  
  #sbatch --export=chr=$i kunkle_fixed_0224_plink_pT.sh
  #sbatch --export=chr=$i plink_wightman_pT.sh
  #sbatch --export=chr=$i kunkle_pT_new_beta.sh
done
