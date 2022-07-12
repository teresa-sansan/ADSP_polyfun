#!/bin/bash
#SBATCH --job-name=pt_with_clump_maxsnp10
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=40G
#SBATCH --time=24:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/cT/chr_sep_plink/wightman/new_beta/%x_%j.log

which_sumstat='wightman'
#clump_path='/gpfs/commons/home/tlin/output/cT/chr_sep_plink/wightman'

qc_sumfile='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/wightman_withbeta_qc.tsv'
sumfile='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/Wightman_2021_hg37_withbeta.tsv'

snp='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/wightman_withbeta.pvalue'
qc_snp='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/wightman_withbeta_qc.pvalue'

if false; then
cat /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/Wightman_2021_hg37_withbeta.tsv | cut -f 1,7 > $snp
cat /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/wightman_withbeta_qc.tsv | cut -f 1,7 > $qc_snp
fi

## try new beta (PIP*beta)
if true; then

clump_path='/gpfs/commons/home/tlin/output/cT/chr_sep_plink/wightman/new_beta'
summary_stat='/gpfs/commons/home/tlin/output/wightman/fixed_0224/finemap/max_snp_10/agg_all_new_beta.tsv'  ## new summary stat.

echo "run chr " $chr
# if [ $chr == 1 ] then;
#   cat $summary_stat |tail -n+2 | cut -f 2 > $clump_path/new_beta/max_snp_10.valid.snp
#   cat agg_all_new_beta.tsv.gz| cut -f 2,16 > ${summary_stat}.pvalue  ## create p value file for that summary stat
# fi

 pvalue=$(echo ${summary_stat}.pvalue)

 awk 'NR!=1{print $3}' ${clump_path}/wightman_${chr}.clumped > ${clump_path}/wightman_${chr}.valid.snp ## only include this line if you did clumping beforehand

## didn't include --extract --q-score-range because I only want to test out the most basic PRS calculation.

 ~/plink \
 --bfile /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/plink_biallelic/filt/${chr}_filt \
 --score  $summary_stat 2 4 16 header \
 --q-score-range range_list.txt $pvalue \
 --extract ${clump_path}/wightman_${chr}.valid.snp \
 --out $clump_path/max_snp_10_clump_pT_chr${chr}
fi


##qc on both
if false; then
#qc='qc'
qc='qc_on_variant_sumstat'
awk 'NR!=1{print $3}' $clump_path/$qc/${which_sumstat}_qc_on_both_clump_chr${chr}.clumped > $clump_path/$qc/chr${chr}.valid.snp
~/plink \
--bfile /gpfs/commons/groups/knowles_lab/data/ADSP_reguloMsL/ADSP_vcf/compact_filtered_vcf_16906/plink_biallelic/qc_on_variant/ADSP_qc_chr${chr} \
--score $qc_sumfile 1 4 9 header \
--q-score-range range_list.txt $qc_snp \
--extract $clump_path/$qc/chr${chr}.valid.snp \
--out $clump_path/$qc/$which_sumstat_pT_chr${chr}
fi

##qc on base
if false; then
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
if false; then
for qc in qc_on_variant qc_on_individual
do
#qc='qc_on_target'
awk 'NR!=1{print $3}' $clump_path/$qc/${which_sumstat}_clump_chr${chr}.clumped > $clump_path/$qc/chr${chr}.valid.snp

~/plink \
--bfile /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/plink_biallelic/$qc/ADSP_qc_chr${chr} \
--score $sumfile 1 4 9 header \
--q-score-range range_list.txt $snp \
--extract $clump_path/$qc/chr${chr}.valid.snp \
--out $clump_path/$qc/$which_sumstat_pT_chr${chr}
done
fi

##no qc
if false; then
qc='before_qc'
awk 'NR!=1{print $3}' $clump_path/$qc/${which_sumstat}_clump_chr${chr}.clumped > $clump_path/$qc/chr${chr}.valid.snp

~/plink \
--bfile /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/plink_biallelic/$qc/ADSP_qc_chr${chr} \
--score $sumfile 1 4 9 header \
--q-score-range range_list.txt $snp \
--extract $clump_path/$qc/chr${chr}.valid.snp \
--out $clump_path/$qc/$which_sumstat_pT_chr${chr}

fi
## SNP_ID, effective allele info, effect size estimate

