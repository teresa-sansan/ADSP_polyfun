#!/bin/bash

#for anno in bl_brain_atac bl_dl_annotations bl
if true; then
for anno in all_except_enformer bl no_ml_new 
do
        for i in 2
        do
                echo "max_snp = $i" 
                sbatch --export=max_snp=$i,anno=$anno  aggregate_finemap_result.sh 
                #sbatch --export=max_snp=$i aggregate_finemap_result.sh 
        done
done
fi


## run converge
if false; then
for anno in bl susie all_anno 
do
sbatch --export=anno=$anno  aggregate_finemap_result.sh 
done
fi