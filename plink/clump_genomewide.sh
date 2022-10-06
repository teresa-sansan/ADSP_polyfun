#!/bin/bash
#SBATCH --job-name=bellenguez_cT
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=70G
#SBATCH --time=21:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/cT/new_plink_genomewide/bellenguez/new_sep22/%x_%j.log


## no qc
plink_path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/plink/'

if false; then
sumstats='wightman'
#sumstats_file='wightman_fixed_beta.tsv'
sumstats_file='wightman_chr_sep/wightman_fixed_beta_chr'
sumstats_qc_file='wightman_chr_sep/wightman_fixed_beta_qc_chr'
output_repo='fixed_rsid_1002'
fi

if true; then
sumstats='bellenguez'
sumstats_file='Bellenguez_et_al_2021_hg37_new_sep20.tsv'
sumstats_qc_file='Bellenguez_et_al_2021_hg37_new_sep20_qc.tsv'
output_repo='new_sep22'
fi

if false; then
for ADSP in ADSP
do
~/plink \
--bfile $plink_path/$ADSP/ADSP_${chr} \
--clump-p1 1 \
--clump-r2 0.1  \
--clump-kb 250  \
--clump /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/${sumstats_file} \
--clump-snp-field SNP \
--clump-field P \
--out /gpfs/commons/home/tlin/output/cT/new_plink_genomewide/$sumstats/$output_repo/$ADSP/$sumstats_${ADSP}_${chr}
done
fi


##qc_all

if true; then
for ADSP in ADSP_qc_all
do
~/plink \
--bfile $plink_path/$ADSP/${ADSP}_${chr} \
--clump-p1 1 \
--clump-r2 0.1  \
--clump-kb 250  \
--clump /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/${sumstats_qc_file} \
--clump-snp-field SNP \
--clump-field P \
--out /gpfs/commons/home/tlin/output/cT/new_plink_genomewide/$sumstats/$output_repo/$ADSP/${ADSP}_${chr} 
done

fi
