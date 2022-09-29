#!/bin/bash
#SBATCH --job-name=pT_wightman
#SBATCH --mail-type=FAIL 
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=30G
#SBATCH --time=40:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/cT/new_plink_genomewide/wightman/ADSP_qc_all/%x_%j.log


plink_path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/plink'


##wightman
if true; then
sumstats='wightman'
beta=10
#bfile='ADSP_qc_all/ADSP_qc_all'
#sumstats_file='wightman_fixed_beta_qc'

bfile='ADSP/ADSP'
sumstats_file='wightman_fixed_beta'
fi

if true; then
#for bfile in ADSP_UKBB_qc/ADSP_UKBB_qc
#for bfile in ADSP_qc_all/ADSP_qc_all ADSP_qc_variant/ADSP_qc_variant ADSP_UKBB/ADSP_UKBB_only ADSP_UKBB_qc/ADSP_UKBB_qc
#for bfile in ADSP/ADSP_all ADSP_qc_all/ADSP_qc_all ADSP_qc_variant/ADSP_qc_variant ADSP_UKBB/ADSP_UKBB_only ADSP_UKBB_qc/ADSP_UKBB_qc


for chr in {1..22}
do

awk 'NR!=1{print $3}' /gpfs/commons/home/tlin/output/cT/new_plink_genomewide/$sumstats/${bfile}_${chr}.clumped > /gpfs/commons/home/tlin/output/cT/new_plink_genomewide/$sumstats/${bfile}_${chr}.valid.snp

~/plink \
--bfile $plink_path/${bfile}_${chr} \
--score  /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/${sumstats_file}.tsv 1 4 $beta header \
--q-score-range range_list.txt /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/${sumstats_file}.snp  \
--extract /gpfs/commons/home/tlin/output/cT/new_plink_genomewide/$sumstats/${bfile}_${chr}.valid.snp \
--out /gpfs/commons/home/tlin/output/cT/new_plink_genomewide/$sumstats/${bfile}_pT_${chr} 

done
fi

