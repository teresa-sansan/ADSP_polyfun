#!/bin/bash
#SBATCH --job-name=36k_clump
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=30G
#SBATCH --time=10:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/prs/pT_36k/wightman/%x_%j.log


## no qc
#plink_path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/plink'

## qc
#plink_path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/plink/ADSP_qc_all'
#bfile='ADSP_no_apoe'


## qc_36k
plink_path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_preview/plink_hg38_qc'
bfile='ADSP.finqc.chr'


# sumstats_file_name='Jansen_qc.tsv'
#sumstats_file_name='Jansen_et_al_2019_hg37_ldsc.tsv'

for chr in {21..22}
do
    ~/plink \
    --bfile /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_preview/plink_hg38_qc/ADSP_qc_chr${chr}\
    --clump-p1 1 \
    --clump-r2 0.1  \
    --clump-kb 250  \
    --clump /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/wightman_fixed_beta_qc.tsv \
    --clump-snp-field SNP \
    --clump-field P \
    --out /gpfs/commons/home/tlin/output/prs/pT_36k/wightman/ADSP_qc_plink_${chr} 

done
