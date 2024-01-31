#!/bin/bash

for i in 3 4 22
do
   echo "run chr $i"
   #sbatch --export=chr=$i /gpfs/commons/home/tlin/script/plink/plink_clump.sh
   #sbatch --export=chr=$i /gpfs/commons/home/tlin/script/plink/wightman_clump.sh
   #sbatch --export=chr=$i /gpfs/commons/home/tlin/script/plink/kunkle_clump_new_beta.sh
   #sbatch --export=chr=$i /gpfs/commons/home/tlin/script/plink/clump_new_beta.sh
   #sbatch --export=chr=$i /gpfs/commons/home/tlin/script/plink/polypred.clump.sh
   #sbatch --export=chr=$i /gpfs/commons/home/tlin/script/plink/kunkle_clump.sh 
   sbatch --export=chr=$i /gpfs/commons/home/tlin/script/plink/clump.sh 

done


