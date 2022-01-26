#!/bin/bash   
#SBATCH --job-name=clump_height
#SBATCH --mail-type=FAIL      
#SBATCH --mail-user=tlin@nygenome.org 
#SBATCH --mem=150G 
#SBATCH --output=/gpfs/commons/home/tlin/data/plink_tutorial/height_test/%x_%j.log


~/plink \
    --bfile /gpfs/commons/home/tlin/data/biallelic/check/${chr}_filt \
    --clump-p1 1 \
    --clump-r2 0.1 \
    --clump-kb 250 \
    --clump /gpfs/commons/home/tlin/data/plink_tutorial/Height.QC.Transformed \
    --clump-snp-field SNP \
    --clump-field P \
    --out /gpfs/commons/home/tlin/data/plink_tutorial/height_test/height_chr${chr}
