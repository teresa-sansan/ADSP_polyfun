#!/bin/bash

for i in {1..22}
do
   echo "run chr $i"
   #sbatch --export=chr=$i /gpfs/commons/home/tlin/script/plink_cT/plink_clump_target.sh
   sbatch --export=chr=$i /gpfs/commons/home/tlin/script/plink_cT/kunkle_clump.sh
done


