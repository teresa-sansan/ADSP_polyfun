#!/bin/bash


for i in  1 5 10
do
	sbatch --export=max_snp=$i polypred.sh
done

