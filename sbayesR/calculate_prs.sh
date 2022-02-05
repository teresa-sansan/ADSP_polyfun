#!/bin/bash
#SBATCH --job-name=sbayesR_prs
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --time=10:00:00
#SBATCH --mem=150G
#SBATCH --output=/gpfs/commons/home/tlin/output/sbayesR%x_%j.log


~/plink \
--bfile /gpfs/commons/home/tlin/data/biallelic/check/${chr}_filt \
--score /gpfs/commons/home/tlin/output/sbayesR/bellenguez_chr${chr}.snpRes 2 5 8 header \
--out /gpfs/commons/home/tlin/output/sbayesR/${chr}_prs
