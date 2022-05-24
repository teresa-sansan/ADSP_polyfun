#!/bin/bash

for chr in {1..4}
do
	sbatch --export=chr=$chr extract_vcf_position.sh
	#sbatch --export=chr=$chr extract_vcf_position38.sh
done
