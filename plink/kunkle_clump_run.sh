#!/bin/bash
#SBATCH --job-name=kunkle_clump
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=150G
#SBATCH --time=10:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/cT/kunkle/fixed_0224/qc_check_target/qc_on_individual/%x_%j.log

##kunkle_QCed
sumstats='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/Kunkle_et_al_2019_hg37_ldsc_qc.tsv.gz'
#sumstats='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/kunkle_2019/Kunkle_etal_Stage1_qc.gz'

## Pvalue, SNP

##qc 


if true; then
echo qc
~/plink \
--bfile /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/plink_biallelic/qc_on_individual/ADSP_qc_chr${chr} \
--clump-p1 1 \
--clump-r2 0.1  \
--clump-kb 250  \
--clump $sumstats \
--clump-snp-field SNP \
--clump-field P \
--out /gpfs/commons/home/tlin/output/cT/kunkle/fixed_0224/qc_check_target/qc_on_individual/kunkle_qc_both_clump_chr${chr}
echo
fi

##qc on base
if false; then
echo qc on base
~/plink \
--bfile ~/data/biallelic/${chr}_filt \
--clump-p1 1 \
--clump-r2 0.1  \
--clump-kb 250  \
--clump $sumstats \
--clump-snp-field SNP \
--clump-field P \
--out /gpfs/commons/home/tlin/output/cT/kunkle/fixed_0224/qc_on_base/kunkle_clump_chr${chr}
echo
fi

## qc on target
sumstats='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/Kunkle_et_al_2019_hg37_ldsc.tsv.gz'
if true; then
echo qc on target
~/plink \
--bfile /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/plink_biallelic/qc_on_individual/ADSP_qc_chr${chr} \
--clump-p1 1 \
--clump-r2 0.1  \
--clump-kb 250  \
--clump $sumstats \
--clump-snp-field SNP \
--clump-field P \
--out /gpfs/commons/home/tlin/output/cT/kunkle/fixed_0224/qc_check_target/qc_on_individual/kunkle_clump_chr${chr}
echo
fi

##before qc
if false; then
echo qc on base
~/plink \
--bfile ~/data/biallelic/${chr}_filt \
--clump-p1 1 \
--clump-r2 0.1  \
--clump-kb 250  \
--clump $sumstats \
--clump-snp-field SNP \
--clump-field P \
--out /gpfs/commons/home/tlin/output/cT/kunkle/fixed_0224/before_qc/kunkle_clump_chr${chr}
fi

