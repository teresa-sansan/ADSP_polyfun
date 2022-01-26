#!/bin/bash
#SBATCH --job-name=height_pT
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=50G
#SBATCH --time=2:00:00
#SBATCH --output=/gpfs/commons/home/tlin/data/plink_tutorial/height_test/only_change_PLINK_input/%x_%j.log

path='/gpfs/commons/home/tlin/data'
clump_path='/gpfs/commons/home/tlin/data/plink_tutorial'

awk 'NR!=1{print $3}' $clump_path/height_test/height_chr${chr}.clumped >  $clump_path/height_test/height_chr${chr}.valid.snp

~/plink \
--bfile /gpfs/commons/home/tlin/data/biallelic/check/${chr}_filt  \
--score ${clump_path}/Height.QC.Transformed 3 4 12 header  \
--q-score-range range_list2.txt  ${clump_path}/SNP.pvalue \
--extract $clump_path/height_test/height_chr${chr}.valid.snp \
--out ${clump_path}/height_test/only_change_PLINK_input/pT_${chr}




##$3 for SNP_ID,$4 for effective allele info, $7 for effect size estimate 

