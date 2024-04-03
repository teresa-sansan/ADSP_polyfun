#!/bin/bash
#SBATCH --job-name=pT_wightman_ibd
#SBATCH --mail-type=FAIL 
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=40G
#SBATCH --time=14:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/prs/pT_36k_ibd/%x_%j.log


#plink_path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/plink'
#plink_path_36k='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_preview/plink_hg38_qc/ADSP_qc_chr'
plink_path_36k='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_QC/annotated_hg37_plink_ibd/qc/qc_chr22'
#clump_path='/gpfs/commons/home/tlin/output/prs/pT_36k/wightman'
clump_path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_QC/annotated_hg37_plink_ibd/plink_prs/clump_wightman'

##wightman
if true; then
sumstats='wightman'
beta=10
#bfile='ADSP_qc_all/ADSP_qc_all'
#sumstats_file='wightman_fixed_beta_qc'

bfile='ADSP/ADSP'
sumstats_file='wightman_fixed_beta'
fi

for chr in {17..22}
do
awk 'NR!=1{print $3}' $clump_path/ADSP_qc_plink_${chr}.clumped > $clump_path/ADSP_qc_plink_${chr}.valid.snp

~/plink \
--bfile $plink_path_36k${chr} \
--score /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/${sumstats_file}.tsv 1 4 $beta header \
--q-score-range range_list.txt /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/${sumstats_file}.snp  \
--extract $clump_path/ADSP_qc_plink_${chr}.valid.snp \
--out $clump_path/ADSP_qc_plink_pT_${chr} 
done

