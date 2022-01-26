#!/bin/bash   
#SBATCH --job-name=clump_height_change_basefile
#SBATCH --mail-type=FAIL      
#SBATCH --mail-user=tlin@nygenome.org 
#SBATCH --mem=150G 
#SBATCH --output=/gpfs/commons/home/tlin/data/plink_tutorial/height_test/change_base_file/%x_%j.log


~/plink \
    --bfile /gpfs/commons/home/tlin/data/plink_tutorial/EUR.QC \
    --clump-p1 1 \
    --clump-r2 0.1 \
    --clump-kb 250 \
    --clump /gpfs/commons/home/tlin/data/bellenguez_2021_final.tsv.gz \
    --clump-snp-field RSID \
    --clump-field PVALUE \
    --out /gpfs/commons/home/tlin/data/plink_tutorial/height_test/change_base_file/height
