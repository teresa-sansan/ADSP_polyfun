#!/bin/bash
for anno in  bl bl_brain_atac bl_dl_annotations susie
do
for i in  1 5 
do
        echo "max_snp = $i" 
        #sbatch --export=max_snp=$i aggregate_finemap_result.sh 
        sbatch --export=max_snp=$i,anno=$anno  aggregate_finemap_result.sh 
	#sbatch --export=max_snp=$i aggregate_finemap_result_kunkle.sh
done
done
