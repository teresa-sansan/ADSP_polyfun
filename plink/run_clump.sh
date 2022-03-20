#!/bin/bash

for i in {1..22}
do
   echo "run chr $i"
   sbatch --export=chr=$i /gpfs/commons/home/tlin/script/plink/plink_clump.sh
   #sbatch --export=chr=$i /gpfs/commons/home/tlin/script/plink/kunkle_clump.sh
done


