#!/bin/bash
#SBATCH --job-name=pT
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=50G
#SBATCH --time=2:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/cT/kunkle/fixed_0224/subsets/%x_%j.log

clump_path='/gpfs/commons/home/tlin/output/cT/kunkle/fixed_0224'

##qc on both
#summary_stat='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/Kunkle_et_al_2019_hg37_ldsc_qc.tsv'
#pvalue='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/Kunkle_et_al_2019_hg37_ldsc_qc.pvalue'

## no apoe
summary_stat='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/Kunkle_remove_APOE_qc.tsv'
pvalue='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/Kunkle_remove_APOE_qc.pvalue'

if true; then
#--bfile /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/plink_biallelic/qc/ADSP_qc_chr${chr} \
#qc='qc'
qc='subsets/qc_on_variant_sumstat'

#awk 'NR!=1{print $3}' $clump_path/$qc/kunkle_qc_both_clump_chr${chr}.clumped > $clump_path/$qc/qc_both_chr${chr}.valid.snp
awk 'NR!=1{print $3}' $clump_path/$qc/kunkle_no_APOE_qc_chr${chr}.clumped > $clump_path/$qc/qc_both_chr${chr}.valid.snp

~/plink \
--bfile /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/plink_biallelic/qc_on_variant/ADSP_qc_chr${chr} \
--score  $summary_stat 1 4 6 header \
--q-score-range range_list.txt $pvalue \
--extract $clump_path/$qc/qc_both_chr${chr}.valid.snp \
--out $clump_path/$qc/qc_variant_sumstat_pT_chr${chr}
fi

##qc on base
if false; then
qc='qc_on_base'
awk 'NR!=1{print $3}' $clump_path/$qc/kunkle_clump_chr${chr}.clumped > $clump_path/$qc/chr${chr}.valid.snp

~/plink \
--bfile /gpfs/commons/home/tlin/data/biallelic/${chr}_filt \
--score  $summary_stat 1 4 6 header \
--q-score-range range_list.txt $pvalue \
--extract $clump_path/$qc/chr${chr}.valid.snp \
--out /gpfs/commons/home/tlin/output/cT/kunkle/fixed_0224/qc_on_base/kunkle_pT_chr${chr}
fi

##qc_on_target
summary_stat='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/Kunkle_et_al_2019_hg37_ldsc.tsv'
pvalue='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/Kunkle_et_al_2019_hg37_ldsc.pvalue'
if false; then
#qc='qc_on_target"

for qc in qc_on_variant qc_on_individual
do
awk 'NR!=1{print $3}' $clump_path/$qc/kunkle_clump_chr${chr}.clumped > $clump_path/$qc/chr${chr}.valid.snp
#--bfile /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/plink_biallelic/qc/ADSP_qc_chr${chr} \

~/plink \
--bfile /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/plink_biallelic/$qc/ADSP_qc_chr${chr} \
--score $summary_stat 1 4 6 header \
--q-score-range range_list.txt $pvalue \
--extract $clump_path/$qc/chr${chr}.valid.snp \
--out /gpfs/commons/home/tlin/output/cT/kunkle/fixed_0224/qc_check_target/$qc/kunkle_pT_chr${chr}
done
fi

##no qc

summary_stat='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/Kunkle_remove_APOE.tsv'
pvalue='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/Kunkle_remove_APOE.pvalue'

if true; then
qc='subsets/before_qc'
#awk 'NR!=1{print $3}' $clump_path/$qc/kunkle_clump_chr${chr}.clumped > $clump_path/$qc/chr${chr}.valid.snp
awk 'NR!=1{print $3}' $clump_path/$qc/kunkle_no_APOE_clump_chr${chr}.clumped > $clump_path/$qc/chr${chr}.valid.snp

~/plink \
--bfile /gpfs/commons/home/tlin/data/biallelic/${chr}_filt \
--score $summary_stat 1 4 6 header \
--q-score-range range_list.txt $pvalue \
--extract $clump_path/$qc/chr${chr}.valid.snp \
--out $clump_path/$qc/kunkle_pT_chr${chr}

fi

##$1 for SNP_ID,$4 for effective allele info, $6 for effect size estimate 

