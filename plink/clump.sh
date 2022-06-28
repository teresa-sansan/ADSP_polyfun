#!/bin/bash
#SBATCH --job-name=clump_genomewide
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
~/plink \
--bfile $plink_path/${bfile}_${chr} \
--clump-p1 1 \
--clump-r2 0.1  \
--clump-kb 250  \
--clump $sumstats_path \
--clump-snp-field SNP \
--clump-field P \
--out /gpfs/commons/home/tlin/output/cT/genomewide_plink/$sumstats/${bfile}_qc_${chr} 

done
fi

