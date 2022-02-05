#!/bin/bash
#SBATCH --job-name=pT
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=50G
#SBATCH --time=2:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/cT/bellenguez/%x_%j.log

path='/gpfs/commons/home/tlin/data'
clump_path='/gpfs/commons/home/tlin/output/cT/bellenguez'
awk 'NR!=1{print $3}' $clump_path/update_RSID_clump_chr${chr}.clumped > $clump_path/clump_chr${chr}.valid.snp
 
~/plink \
--bfile /gpfs/commons/home/tlin/data/biallelic/check/${chr}_filt \
--score /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/bellenguez_2021/bellenguez_2021_final_rename.tsv 3 4 7 header \
--q-score-range range_list.txt /gpfs/commons/home/tlin/data/bellenguez_snp.pvalue \
--extract $clump_path/clump_chr${chr}.valid.snp \
--out /gpfs/commons/home/tlin/output/cT/bellenguez/bellenguez_pT_chr${chr}



##kunkle
#clump_path='/gpfs/commons/home/tlin/output/cT/kunkle/'
#clump_path='/gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap/clumping'
#--score ${path}/kunkle_etal_Stage1_info.tsv 3 4 8 header \
#--q-score-range range_list.txt $path/kunkle_check.pvalue \
#--out /gpfs/commons/home/tlin/output/cT/kunkle/kunkle_pT_chr${chr}


##$3 for SNP_ID,$4 for effective allele info, $7 for effect size estimate 

