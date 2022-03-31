#!/bin/bash
#SBATCH --job-name=pT
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=50G
#SBATCH --time=2:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/cT/wightman/%x_%j.log

which_sumstat='wightman'
clump_path='/gpfs/commons/home/tlin/output/cT/wightman'
qc_sumfile='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/wightman_withbeta_qc.tsv'
sumfile='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/Wightman_2021_hg37_withbeta.tsv'

snp='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/wightman_withbeta.pvalue'
qc_snp='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/wightman_withbeta_qc.pvalue'

if false; then
cat /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/Wightman_2021_hg37_withbeta.tsv | cut -f 1,7 > $snp
cat /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/wightman_withbeta_qc.tsv | cut -f 1,7 > $qc_snp
fi

##qc on both
if true; then
qc='qc'

awk 'NR!=1{print $3}' $clump_path/$qc/${which_sumstat}_clump_chr${chr}.clumped > $clump_path/$qc/chr${chr}.valid.snp
~/plink \
--bfile /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/plink_biallelic/qc/ADSP_qc_chr${chr} \
--score $qc_sumfile 1 4 9 header \
--q-score-range range_list.txt $qc_snp \
--extract $clump_path/$qc/chr${chr}.valid.snp \
--out $clump_path/$qc/$which_sumstat_pT_chr${chr}
fi

##qc on base
if true; then
qc='qc_on_base'
awk 'NR!=1{print $3}' $clump_path/$qc/${which_sumstat}_clump_chr${chr}.clumped > $clump_path/$qc/chr${chr}.valid.snp
~/plink \
--bfile /gpfs/commons/home/tlin/data/biallelic/${chr}_filt \
--score $qc_sumfile 1 4 9 header \
--q-score-range range_list.txt $qc_snp \
--extract $clump_path/$qc/chr${chr}.valid.snp \
--out $clump_path/$qc/$which_sumstat_pT_chr${chr}
fi

##qc_on_target
if true; then
qc='qc_on_target'
awk 'NR!=1{print $3}' $clump_path/$qc/${which_sumstat}_clump_chr${chr}.clumped > $clump_path/$qc/chr${chr}.valid.snp

~/plink \
--bfile /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/plink_biallelic/qc/ADSP_qc_chr${chr} \
--score $sumfile 1 4 9 header \
--q-score-range range_list.txt $snp \
--extract $clump_path/$qc/chr${chr}.valid.snp \
--out $clump_path/$qc/$which_sumstat_pT_chr${chr}
fi

##no qc
if true; then
qc='before_qc'
awk 'NR!=1{print $3}' $clump_path/$qc/${which_sumstat}_clump_chr${chr}.clumped > $clump_path/$qc/chr${chr}.valid.snp

~/plink \
--bfile /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/plink_biallelic/qc/ADSP_qc_chr${chr} \
--score $sumfile 1 4 9 header \
--q-score-range range_list.txt $snp \
--extract $clump_path/$qc/chr${chr}.valid.snp \
--out $clump_path/$qc/$which_sumstat_pT_chr${chr}

fi
##$1 for SNP_ID,$4 for effective allele info, $9 for effect size estimate 

