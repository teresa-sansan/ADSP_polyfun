#!/bin/bash
#SBATCH --job-name=jansen_clump
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=30G
#SBATCH --time=10:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/cT/new_plink_genomewide/jansen/ADSP_qc_all/%x_%j.log


## no qc
#plink_path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/plink'

## qc
plink_path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/plink/ADSP_qc_all'
bfile='ADSP_no_apoe'

if true; then
#for bfile in ADSP_UKBB_qc/ADSP_UKBB_qc
#for bfile in ADSP_qc_all/ADSP_qc_all ADSP_qc_variant/ADSP_qc_variant ADSP_UKBB/ADSP_UKBB_only ADSP_UKBB_qc/ADSP_UKBB_qc
#for bfile in ADSP/ADSP_all ADSP_qc_all/ADSP_qc_all ADSP_qc_variant/ADSP_qc_variant ADSP_UKBB/ADSP_UKBB_only ADSP_UKBB_qc/ADSP_UKBB_qc


sumstats_file_name='Jansen_qc.tsv'
#sumstats_file_name='Jansen_et_al_2019_hg37_ldsc.tsv'

~/plink \
--bfile $plink_path/ADSP_qc_all_${chr} \
--clump-p1 1 \
--clump-r2 0.1  \
--clump-kb 250  \
--clump /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/${sumstats_file_name} \
--clump-snp-field SNP \
--clump-field P \
--out /gpfs/commons/home/tlin/output/cT/new_plink_genomewide/jansen/ADSP_qc_all/ADSP_qc_plink_${chr} 

fi

