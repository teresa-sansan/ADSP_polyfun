#!/bin/bash
#SBATCH --job-name=wightman_clump_qc_variant_sumstat
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=30G
#SBATCH --time=10:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/cT/wightman/qc_on_variant_sumstat/%x_%j.log

##wightman
sumstats='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/wightman_qc.tsv.gz'
output='/gpfs/commons/home/tlin/output/cT/wightman'
## Pvalue, SNP

##qc
#--bfile /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/plink_biallelic/qc/ADSP_qc_chr${chr} \

if true; then
echo start chr $chr
~/plink \
--bfile /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/plink_biallelic/qc_on_variant/ADSP_qc_chr${chr} \
--clump-p1 1 \
--clump-r2 0.1  \
--clump-kb 250  \
--clump $sumstats \
--clump-snp-field SNP \
--clump-field P \
--out /gpfs/commons/home/tlin/output/cT/wightman/qc_on_variant_sumstat/wightman_qc_on_both_clump_chr${chr}

fi

#--out /gpfs/commons/home/tlin/output/cT/wightman/qc/wightman_clump_chr${chr}

##qc_base
if false; then
echo start chr $chr
~/plink \
--bfile ~/data/biallelic/${chr}_filt \
--clump-p1 1 \
--clump-r2 0.1  \
--clump-kb 250  \
--clump $sumstats \
--clump-snp-field SNP \
--clump-field P \
--out /gpfs/commons/home/tlin/output/cT/wightman/qc_on_base/wightman_clump_chr${chr}
fi

sumstats='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/Bellenguez_et_al_2021_hg37_no_dup.tsv.gz'
#sumstats='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/Wightman_et_al_2021_hg37_ldsc.tsv.gz'


##qc_target
if false; then
for qc in qc_on_variant qc_on_individual
do
echo start chr $chr
~/plink \
--bfile /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/plink_biallelic/$qc/ADSP_qc_chr${chr} \
--clump-p1 1 \
--clump-r2 0.1  \
--clump-kb 250  \
--clump $sumstats \
--clump-snp-field SNP \
--clump-field $pvalue \
--out $output/$qc/bellenguez_clump_chr${chr}
done

fi

##before qc
if false; then
echo start chr $chr
~/plink \
--bfile ~/data/biallelic/${chr}_filt \
--clump-p1 1 \
--clump-r2 0.1  \
--clump-kb 250  \
--clump $sumstats \
--clump-snp-field SNP \
--clump-field P \
--out /gpfs/commons/home/tlin/output/cT/wightman/before_qc/wightman_clump_chr${chr}
fi
