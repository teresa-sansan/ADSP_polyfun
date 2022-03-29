#!/bin/bash

for i in {1..22}
do
  echo "start calculating PRS for chr $i"
  #sbatch --export=chr=$i plink_pT.sh   
  #sbatch --export=chr=$i kunkle_plink_pT.sh   
  #sbatch --export=chr=$i kunkle_fixed_0224_plink_pT.sh
  sbatch --export=chr=$i plink_wightman_pT.sh
done
