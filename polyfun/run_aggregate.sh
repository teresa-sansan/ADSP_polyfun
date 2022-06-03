#!/bin/bash

for i in  1 3 5 7 10
do
        echo "max_snp = $i" 
        #sbatch --export=max_snp=$i aggregate_finemap_result.sh 
	sbatch --export=max_snp=$i aggregate_finemap_result_kunkle.sh
done

