#!/bin/bash
#SBATCH --job-name=kunkle_clump
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=150G
#SBATCH --time=10:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/cT/new_plink_genomewide/kunkle/ADSP_no_apoe/%x_%j.log

##kunkle_QCed
#sumstats='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/Kunkle_et_al_2019_hg37_ldsc_qc.tsv.gz'
#sumstats='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/kunkle_2019/Kunkle_etal_Stage1_qc.gz'


## kunkle_no_APOE
sumstats='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/Kunkle_remove_APOE_qc.tsv'


## kunkle_polyfun
#sumstats='/gpfs/commons/home/tlin/output/kunkle/kunkle_fixed_0224/finemap/max_snp_10/aggregate.all.txt.gz'

## Pvalue, SNP
##qc on both
#--bfile /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/plink_biallelic/qc_on_variant/ADSP_qc_chr${chr} \

if true; then
echo qc
~/plink \
--bfile /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/plink/ADSP_qc_all/ADSP_qc_all_${chr} \
--clump-p1 1 \
--clump-r2 0.1  \
--clump-kb 250  \
--clump $sumstats \
--clump-snp-field SNP \
--clump-field P \
--out /gpfs/commons/home/tlin/output/cT/new_plink_genomewide/kunkle/ADSP_no_apoe/kunkle_noAPOE_qc_chr${chr}
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
sumstats='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/Kunkle_et_al_2019_hg37_ldsc.tsv'
if false; then
for qc in  qc_on_variant qc_on_individual
do
~/plink \
--bfile /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/plink_biallelic/$qc/ADSP_qc_chr${chr} \
--clump-p1 1 \
--clump-r2 0.1  \
--clump-kb 250  \
--clump $sumstats \
--clump-snp-field SNP \
--clump-field P \
--out /gpfs/commons/home/tlin/output/cT/kunkle/fixed_0224/qc_check_target/$qc/kunkle_clump_chr${chr}
done
fi

sumstats='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/Kunkle_remove_APOE_qc.tsv'
##before qc
if false; then
echo no qc
~/plink \
--bfile ~/data/biallelic/${chr}_filt \
--clump-p1 1 \
--clump-r2 0.1  \
--clump-kb 250  \
--clump $sumstats \
--clump-snp-field SNP \
--clump-field P \
--out /gpfs/commons/home/tlin/output/cT/kunkle/fixed_0224/subsets/before_qc/kunkle_no_APOE_clump_chr${chr}
fi
