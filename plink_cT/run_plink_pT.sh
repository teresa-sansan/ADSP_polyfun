#!/bin/bash


for i in {1..22}
do
  echo "start calculating PRS for chr $i"
  sbatch --export=chr=$i plink_pT.sh   
  #sbatch --export=chr=$i /gpfs/commons/home/tlin/polyfun_script/test/pT_height.sh
done
