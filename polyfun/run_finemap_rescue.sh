#!/bin/bash 

for anno in bl bl_brain_atac bl_dl_annotations
do
    sbatch --export=max_snp=5,anno=${anno} finemap_rescue.sh 
    sbatch --export=max_snp=10,anno=${anno} finemap_rescue.sh 
done




