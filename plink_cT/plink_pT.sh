#!/bin/bash
#SBATCH --job-name=pT
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=50G
#SBATCH --time=2:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/cT/kunkle/%x_%j.log

path='/gpfs/commons/home/tlin/data'
#clump_path='/gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap/clumping'
clump_path='/gpfs/commons/home/tlin/output/cT/kunkle/'
#awk 'NR!=1{print $3}' $clump_path/clump_chr${chr}.clumped >  $clump_path/clump_chr${chr}.valid.snp
 
~/plink \
--bfile /gpfs/commons/home/tlin/data/biallelic/check/${chr}_filt \
--score ${path}/kunkle_etal_Stage1_info.tsv 3 4 8 header \
--q-score-range range_list.txt $path/kunkle_check.pvalue \
--extract $clump_path/clump_chr${chr}.valid.snp \
--out /gpfs/commons/home/tlin/output/cT/kunkle/kunkle_pT_chr${chr}


#--score ${path}/bellenguez_2021_final.tsv 3 4 7 header \


##$3 for SNP_ID,$4 for effective allele info, $7 for effect size estimate 

