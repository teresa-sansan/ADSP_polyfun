#!/bin/bash
#SBATCH --job-name=pT_wightman
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=50G
#SBATCH --time=10:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/cT/new_plink_genomewide/wightman/fixed_rsid_1002/ADSP_qc_all/%x_%j.log

#clump_path='/gpfs/commons/home/tlin/output/cT/bellenguez/fixed_0224'
#clump_path='/gpfs/commons/home/tlin/output/cT/new_plink_genomewide/kunkle/'

##1 4 6
#sumstat_file='Kunkle_remove_APOE_qc.tsv'

## 1 4 13
#sumstat_file='Jansen_qc.tsv'
#clump_path='/gpfs/commons/home/tlin/output/cT/new_plink_genomewide/jansen/'

## 1 4 10 new wightman
clump_path='/gpfs/commons/home/tlin/output/cT/new_plink_genomewide/wightman/fixed_rsid_1002/'
sumstat_file='wightman_chr_sep/wightman_fixed_beta_qc_chr'
processed='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed'

for chr in {5..22}
do
##qc on both
if true; then
#qc='qc'
#qc='qc_on_variant_sumstat'
#qc='ADSP_no_apoe'
qc='ADSP_qc_all'
awk 'NR!=1{print $3}' $clump_path/$qc/ADSP_qc_all_${chr}.clumped > $clump_path/$qc/chr${chr}.valid.snp
~/plink \
--bfile  /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/plink/ADSP_qc_all/ADSP_qc_all_${chr} \
--score $processed/${sumstat_file}${chr}.tsv 1 4 10 header \
--q-score-range range_list.txt $processed/${sumstat_file}${chr}.pvalue \
--extract $clump_path/$qc/chr${chr}.valid.snp \
--out $clump_path/$qc/${qc}_pT_chr${chr}
fi

#no qc
if false; then
sumstat_file='wightman_chr_sep/wightman_fixed_beta_chr'
qc='ADSP'
awk 'NR!=1{print $3}' $clump_path/$qc/ADSP_${chr}.clumped > $clump_path/$qc/chr${chr}.valid.snp

~/plink \
--bfile /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/plink/ADSP/ADSP_${chr} \
--score $processed/${sumstat_file}${chr}.tsv 1 4 10 header \
--q-score-range range_list.txt $processed/${sumstat_file}${chr}.pvalue  \
--extract $clump_path/$qc/chr${chr}.valid.snp \
--out $clump_path/$qc/${qc}_pT_chr${chr}
fi

##qc on base
if false; then
qc='qc_on_base'
#awk 'NR!=1{print $3}' $clump_path/$qc/bellenguez_clump_chr${chr}.clumped > $clump_path/$qc/chr${chr}.valid.snp

~/plink \
--bfile /gpfs/commons/home/tlin/data/biallelic/${chr}_filt \
--score /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/Bellenguez_et_al_2021_hg37_qc.tsv 1 4 7 header \
--q-score-range range_list.txt /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/Bellenguez_et_al_2021_hg37_no_dup.pvalue \
--extract $clump_path/$qc/chr${chr}.valid.snp \
--out /gpfs/commons/home/tlin/output/cT/bellenguez/fixed_0224/qc_on_base/bellenguez_pT_chr${chr}
fi

##qc_on_target
if false; then
#qc='qc_on_target'
for qc in qc_on_individual qc_on_variant
do
awk 'NR!=1{print $3}' $clump_path/$qc/bellenguez_clump_chr${chr}.clumped > $clump_path/$qc/chr${chr}.valid.snp

~/plink \
--bfile /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/plink_biallelic/$qc/ADSP_qc_chr${chr} \
--score /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/Bellenguez_et_al_2021_hg37_no_dup.tsv  1 4 7 header \
--q-score-range range_list.txt /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/Bellenguez_et_al_2021_hg37_qc.pvalue \
--extract $clump_path/$qc/chr${chr}.valid.snp \
--out /gpfs/commons/home/tlin/output/cT/bellenguez/fixed_0224/$qc/bellenguez_pT_chr${chr}
done
fi

done

#
##$1 for SNP_ID,$4 for effective allele info, $7 for effect size estimate 

