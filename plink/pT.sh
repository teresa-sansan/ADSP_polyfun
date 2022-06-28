#!/bin/bash
#SBATCH --job-name=pT_genomewide_bychr
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=30G
#SBATCH --time=40:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/cT/genomewide_plink/%x_%j.log


## no qc
plink_path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/plink'

if true; then
#for bfile in ADSP_UKBB_qc/ADSP_UKBB_qc
#for bfile in ADSP_qc_all/ADSP_qc_all ADSP_qc_variant/ADSP_qc_variant ADSP_UKBB/ADSP_UKBB_only ADSP_UKBB_qc/ADSP_UKBB_qc
for bfile in ADSP/ADSP_all ADSP_qc_all/ADSP_qc_all ADSP_qc_variant/ADSP_qc_variant ADSP_UKBB/ADSP_UKBB_only ADSP_UKBB_qc/ADSP_UKBB_qc

do

awk 'NR!=1{print $3}' /gpfs/commons/home/tlin/output/cT/genomewide_plink/$sumstats/${bfile}_qc_${chr}.clumped > /gpfs/commons/home/tlin/output/cT/genomewide_plink/$sumstats/${bfile}_qc_${chr}.valid.snp

~/plink \
--bfile $plink_path/${bfile}_${chr} \
--score  $sumstats_path 1 4 $beta header \
--q-score-range range_list.txt /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/$pvalue  \
--extract /gpfs/commons/home/tlin/output/cT/genomewide_plink/$sumstats/${bfile}_qc_${chr}.valid.snp \
--out /gpfs/commons/home/tlin/output/cT/genomewide_plink/$sumstats/${bfile}_pT_${chr} 

done
fi

