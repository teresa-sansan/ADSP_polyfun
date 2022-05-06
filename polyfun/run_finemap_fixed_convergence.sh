#!/bin/bash

for max_snp in 3
do
sbatch --export=max_num_snp=$max_snp finemap_fixed_convergence.sh		
done

