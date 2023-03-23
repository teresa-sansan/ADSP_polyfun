#!/bin/bash 
for i in {1..10}
do
echo run chr $i
sbatch --export=chr=$i maf_biallelic_filter.sh
done