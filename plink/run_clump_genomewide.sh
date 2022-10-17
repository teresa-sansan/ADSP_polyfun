#!/bin/bash

for i in {1..22}
do
sbatch --export=chr=$i /gpfs/commons/home/tlin/script/plink/clump_genomewide.sh

done


