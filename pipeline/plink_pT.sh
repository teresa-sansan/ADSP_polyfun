#!/bin/bash
#SBATCH --job-name=bellenguez_pT
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=50G
#SBATCH --time=2:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap/pT/%x_%j.log

path='/gpfs/commons/home/tlin/data'
 
~/plink \
--bfile /gpfs/commons/home/tlin/data/biallelic/${chr}_filt \
--score ${path}/bellenguez_2021_final.tsv 3 4 7 'header' \
--q-score-range range_list.txt  $path/bellenguez.pvalue \
--extract /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap/clumping/clump_chr${chr}.clumped \
--out /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap/pT/pT_chr${chr}




##$3 for SNP_ID,$4 for effective allele info, $7 for effect size estimate 

