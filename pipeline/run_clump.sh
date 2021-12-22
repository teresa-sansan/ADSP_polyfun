#!/bin/bash

for i in {1..22}
do
   echo "run chr $i"
   sbatch --export=chr=$i plink_clump.sh
done
