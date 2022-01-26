#!/bin/bash
#SBATCH --job-name=height_pT
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=50G
#SBATCH --time=2:00:00
#SBATCH --output=/gpfs/commons/home/tlin/data/plink_tutorial/height_test/change_base_file/%x_%j.log

path='/gpfs/commons/home/tlin/data'
clump_path='/gpfs/commons/home/tlin/data/plink_tutorial/height_test/change_base_file'

awk 'NR!=1{print $3}' $clump_path/height.clumped >  $clump_path/height.valid.snp

~/plink \
--bfile /gpfs/commons/home/tlin/data/plink_tutorial/EUR.QC  \
--score ${path}/bellenguez_2021_final.tsv 3 4 7 header  \
--q-score-range /gpfs/commons/home/tlin/polyfun_script/pipeline/range_list.txt $path/bellenguez.pvalue \
--extract $clump_path/height.valid.snp  \
--out ${clump_path}/height_pT




##$3 for SNP_ID,$4 for effective allele info, $7 for effect size estimate 

