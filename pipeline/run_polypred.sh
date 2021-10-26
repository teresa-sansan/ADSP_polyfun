#!/bin/bash


for i in 1 3 5 7 10
do
	sbatch --export=max_snp=$i polypred.sh
done

