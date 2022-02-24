#!/bin/bash
#SBATCH --job-name=pT
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=50G
#SBATCH --time=2:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/cT/kunkle/qc_check/%x_%j.log


##qc_on_base
clump_path='/gpfs/commons/home/tlin/output/cT/kunkle/qc_on_base'
#awk 'NR!=1{print $3}' $clump_path/kunkle_clump_chr${chr}.clumped > $clump_path/chr${chr}.valid.snp

if false; then
echo run_qc_on_summary_stats
#~/plink \
#--bfile /gpfs/commons/home/tlin/data/biallelic/${chr}_filt \
#--score /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/kunkle_2019/Kunkle_etal_Stage1_qc.tsv 3 4 6 header \
#--q-score-range range_list.txt /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/kunkle_2019/Kunkle_SNP.pvalue \
#--extract $clump_path/chr${chr}.valid.snp \
#--out /gpfs/commons/home/tlin/output/cT/kunkle/qc_on_base/kunkle_pT_chr${chr}
fi



##qc_on_target
clump_path='/gpfs/commons/home/tlin/output/cT/kunkle/qc_on_target'
#awk 'NR!=1{print $3}' $clump_path/kunkle_clump_chr${chr}.clumped > $clump_path/chr${chr}.valid.snp
if false; then
echo run qc_on_target
#~/plink \
#--bfile /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/plink_biallelic/qc/ADSP_qc_chr${chr} \
#--score /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/kunkle_2019/kunkle_etal_Stage1_info.tsv 3 4 8 header \
#--q-score-range range_list.txt /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/kunkle_2019/Kunkle_no_qc.SNP \
#--extract $clump_path/chr${chr}.valid.snp \
#--out /gpfs/commons/home/tlin/output/cT/kunkle/qc_on_target/kunkle_pT_chr${chr}
fi


##qc_on_everything
clump_path='/gpfs/commons/home/tlin/output/cT/kunkle/qc_check'
awk 'NR!=1{print $3}' $clump_path/kunkle_clump_chr${chr}.clumped > $clump_path/chr${chr}.valid.snp
if true; then
echo run_qc_on_everything
~/plink \
--bfile /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/plink_biallelic/qc/ADSP_qc_chr${chr} \
--score /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/kunkle_2019/Kunkle_etal_Stage1_qc.tsv 3 4 6 header \
--q-score-range range_list.txt /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/kunkle_2019/Kunkle_SNP.pvalue \
--extract $clump_path/chr${chr}.valid.snp \
--out /gpfs/commons/home/tlin/output/cT/kunkle/qc_check/kunkle_pT_chr${chr}

fi
 
