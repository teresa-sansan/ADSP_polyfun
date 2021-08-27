#!/bin/bash


for i in {1..22}
do
	echo "run chr$i" 
	sbatch --export=chr=$i baseline_finemap5.sh
done
