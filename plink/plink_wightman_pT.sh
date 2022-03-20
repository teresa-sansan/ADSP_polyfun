#!/bin/bash
#SBATCH --job-name=pT
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=50G
#SBATCH --time=2:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/cT/wightman/%x_%j.log

clump_path='/gpfs/commons/home/tlin/output/cT/wightman'

##qc on both
if true; then
qc='qc'
#awk 'NR!=1{print $3}' $clump_path/$qc/bellenguez_clump_chr${chr}.clumped > $clump_path/$qc/chr${chr}.valid.snp
~/plink \
--bfile /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/plink_biallelic/qc/ADSP_qc_chr${chr} \
--score /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/Bellenguez_et_al_2021_hg37_qc.tsv 1 4 7 header \
--q-score-range range_list.txt /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/Bellenguez_et_al_2021_hg37_qc.pvalue \
--extract $clump_path/$qc/chr${chr}.valid.snp \
--out /gpfs/commons/home/tlin/output/cT/bellenguez/fixed_0224/qc/bellenguez_pT_chr${chr}
fi

##qc on base
if false; then
qc='qc_on_base'
#awk 'NR!=1{print $3}' $clump_path/$qc/bellenguez_clump_chr${chr}.clumped > $clump_path/$qc/chr${chr}.valid.snp

~/plink \
--bfile /gpfs/commons/home/tlin/data/biallelic/${chr}_filt \
--score /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/Bellenguez_et_al_2021_hg37_qc.tsv 1 4 7 header \
--q-score-range range_list.txt /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/Bellenguez_et_al_2021_hg37_no_dup.pvalue \
--extract $clump_path/$qc/chr${chr}.valid.snp \
--out /gpfs/commons/home/tlin/output/cT/bellenguez/fixed_0224/qc_on_base/bellenguez_pT_chr${chr}
fi

##qc_on_target
if false; then
qc='qc_on_target'
#awk 'NR!=1{print $3}' $clump_path/$qc/bellenguez_clump_chr${chr}.clumped > $clump_path/$qc/chr${chr}.valid.snp

~/plink \
--bfile /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/plink_biallelic/qc/ADSP_qc_chr${chr} \
--score /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/Bellenguez_et_al_2021_hg37_no_dup.tsv 1 4 7 header \
--q-score-range range_list.txt /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/Bellenguez_et_al_2021_hg37_no_dup.pvalue \
--extract $clump_path/$qc/chr${chr}.valid.snp \
--out /gpfs/commons/home/tlin/output/cT/bellenguez/fixed_0224/qc_on_target/bellenguez_pT_chr${chr}
fi

##no qc
if false; then
qc='before_qc'
#awk 'NR!=1{print $3}' $clump_path/$qc/bellenguez_clump_chr${chr}.clumped > $clump_path/$qc/chr${chr}.valid.snp

~/plink \
--bfile /gpfs/commons/home/tlin/data/biallelic/${chr}_filt \
--score /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/Bellenguez_et_al_2021_hg37_no_dup.tsv 1 4 7 header \
--q-score-range range_list.txt /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/Bellenguez_et_al_2021_hg37_no_dup.pvalue \
--extract $clump_path/$qc/chr${chr}.valid.snp \
--out /gpfs/commons/home/tlin/output/cT/bellenguez/fixed_0224/before_qc/bellenguez_pT_chr${chr}

fi
##$1 for SNP_ID,$4 for effective allele info, $7 for effect size estimate 

