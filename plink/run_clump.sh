#!/bin/bash

for i in {1..21}
do
   echo "run chr $i"
   #sbatch --export=chr=$i /gpfs/commons/home/tlin/script/plink/plink_clump.sh
   #sbatch --export=chr=$i /gpfs/commons/home/tlin/script/plink/wightman_clump.sh
   #sbatch --export=chr=$i /gpfs/commons/home/tlin/script/plink/kunkle_clump_new_beta.sh
   #sbatch --export=chr=$i /gpfs/commons/home/tlin/script/plink/clump_new_beta.sh
   sbatch --export=chr=$i /gpfs/commons/home/tlin/script/plink/polypred.clump.sh
done


