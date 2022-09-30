#!/bin/bash
#SBATCH --job-name=wightman_new
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=50G
#SBATCH --time=11:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/cT/new_plink_genomewide/wightman/ADSP/%x_%j.log


## no qc
plink_path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/plink/'
sumstats='wightman'
#sumstats_file='wightman_fixed_beta.tsv'
## to test out chr2, try chr only sumstat
sumstats_file='wightman_fixed_beta_chr2.tsv'
sumstats_qc_file='wightman_fixed_beta_qc.tsv'

chr=2

if true; then
for ADSP in ADSP
do
~/plink \
--bfile $plink_path/$ADSP/ADSP_${chr} \
--clump-p1 1 \
--clump-r2 0.1  \
--clump-kb 250  \
--clump /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/$sumstats_file \
--clump-snp-field SNP \
--clump-field P \
--out /gpfs/commons/home/tlin/output/cT/new_plink_genomewide/$sumstats/$ADSP/$sumstats_${ADSP}_${chr}

done
fi


##qc_all

if false; then
for ADSP in ADSP_qc_all
do
~/plink \
--bfile $plink_path/$ADSP/${ADSP}_${chr} \
--clump-p1 1 \
--clump-r2 0.1  \
--clump-kb 250  \
--clump /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/$sumstats_qc_file \
--clump-snp-field SNP \
--clump-field P \
--out /gpfs/commons/home/tlin/output/cT/new_plink_genomewide/$sumstats/$ADSP/${ADSP}_${chr} 

done

fi
