#!/bin/bash
#SBATCH --job-name=wightman_clump
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=5G
#SBATCH --time=25:00:00
#SBATCH --output=/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_QC/annotated_hg37_plink_ibd/plink_prs/clump_bellenguez/%x_%j.log

## qc_36k
plink_path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_QC/annotated_hg37_plink_ibd/qc'
bfile='qc_chr'


wightman='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/wightman_fixed_beta_qc.tsv'
bellenguez='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/Bellenguez_et_al_2021_hg37_new_sep20_qc.tsv'

~/plink \
--bfile $plink_path/$bfile${chr} \
--clump-p1 1 \
--clump-r2 0.1  \
--clump-kb 250  \
--clump $wightman \
--clump-snp-field SNP \
--clump-field P \
--out /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_QC/annotated_hg37_plink_ibd/plink_prs/clump_wightman/ADSP_qc_plink_${chr} 



