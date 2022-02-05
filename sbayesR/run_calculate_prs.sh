#!/bin/bash

for i in {1..22}
do
  echo start running chr $i
  sbatch --export=chr=$i /gpfs/commons/home/tlin/script/sbayesR/calculate_prs.sh

done
