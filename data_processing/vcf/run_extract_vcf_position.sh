#!/bin/bash

for chr in {7..14}
do
	sbatch --export=chr=$chr extract_vcf_position.sh

done
