#!/bin/bash
#SBATCH --job-name=36k_bellenguez_clump_pt_
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=50G
#SBATCH --time=30:00:00
#SBATCH --array=1-22%13
#SBATCH --output=/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/PRS/36k_hg38/clumping_pt_oct/prs_middlefile/%x_%j.log

## qc_36k
plink_path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_QC/annotated_hg37_plink_ibd/qc/qc_chr'
wightman='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/wightman_fixed_beta_qc.tsv'

## 36k_hg38
sumstat_path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/'
#bellenguez='bellenguez_hg38_from_parquet_allinfo.tsv'
bellenguez='bellenguez_hg38_parquet_flipped_a1a2.tsv'
plink_36k_hg38='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_QC/annotated_hg38_plink_qc/ADSP.chr' ## this is for prs

#plink_path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/LD/plink_file_hg38/ADSP_EUR_chr' ## here only use EUR
plink_path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/LD/plinkfile_hg38_rerun/ADSP_EUR_chr'
clump_path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/PRS/36k_hg38/clumping_pt_oct/clump_bellenguez/eur_plink_'
output='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/PRS/36k_hg38/clumping_pt_oct/prs_middlefile/'

chr=$SLURM_ARRAY_TASK_ID

# clumping
# ~/plink \
# --bfile ${plink_path}${chr} \
# --clump-p1 1 \
# --clump-r2 0.5  \
# --clump-kb 250  \
# --clump $sumstat_path/$bellenguez \
# --clump-snp-field SNP \
# --clump-field P \
# --out ${clump_path}${chr} 

# ## get valid snp
# awk 'NR!=1{print $3}' ${clump_path}${chr}.clumped > ${clump_path}${chr}.valid.snp

## p value thresholding
~/plink \
--allow-no-sex \
--bfile ${plink_36k_hg38}${chr} \
--score ${sumstat_path}/${bellenguez} 3 4 6 header \
--out $output/prs_plink_pT_${chr} 

# ## SNP, Effective Allele, BETA
# --q-score-range range_list.txt $sumstat_path/bellenguez_hg38_parquet_flipped_a1a2.snp \
# --extract ${clump_path}${chr}.valid.snp \