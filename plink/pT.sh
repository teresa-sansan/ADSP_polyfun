#!/bin/bash
#SBATCH --job-name=bellenguez_ibd
#SBATCH --mail-type=FAIL 
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=40G
#SBATCH --time=14:00:00
#SBATCH --array=1-22%12
#SBATCH --output=/gpfs/commons/home/tlin/output/prs/pT_36k_ibd/%x_%j.log


#plink_path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/plink'
#plink_path_36k='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_preview/plink_hg38_qc/ADSP_qc_chr'
#plink_path_36k='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_QC/annotated_hg37_plink_ibd/qc/qc_chr'
#clump_path='/gpfs/commons/home/tlin/output/prs/pT_36k/wightman'
#clump_path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_QC/annotated_hg37_plink_ibd/plink_prs/clump_bellenguez'


plink_36k_hg38='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_QC/annotated_hg38_plink_qc/ADSP.chr'
clump_path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/PRS/36k_hg38/c+pt/clump_bellenguez/'



if true; then
sumstats='bellenguez'
beta=11 ## wightman is 10
#bfile='ADSP_qc_all/ADSP_qc_all'
#sumstats_file='wightman_fixed_beta_qc'

bfile='ADSP/ADSP'
sumstats_file='wightman_fixed_beta'
fi

#chr=$SLURM_ARRAY_TASK_ID
chr=$1
awk 'NR!=1{print $3}' $clump_path/ADSP_qc_plink_${chr}.clumped > $clump_path/ADSP_qc_plink_${chr}.valid.snp

~/plink \
--bfile $plink_path_36k${chr} \
--score  /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/Bellenguez_et_al_2021_hg37_new_sep20_qc_nodup.tsv 1 4 $beta header \
--q-score-range range_list_polyfun.txt /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/${sumstats_file}.snp  \
--extract $clump_path/ADSP_qc_plink_${chr}.valid.snp \
--out $clump_path/ADSP_qc_plink_pT_${chr} 

#/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/${sumstats_file}.tsv
#/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/sumstat_chirag 