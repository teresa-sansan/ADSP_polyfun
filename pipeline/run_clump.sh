#!/bin/bash

for i in {1..22}
do
   echo "run chr $i"
   sbatch --export=chr=$i /gpfs/commons/home/tlin/polyfun_script/test/clumping_height.sh
done


   #sbatch --export=chr=$i plink_clump.sh
