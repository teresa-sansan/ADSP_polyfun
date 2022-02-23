#!/bin/bash
#SBATCH --job-name=pT
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=50G
#SBATCH --time=2:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/cT/bellenguez/qc_on_base/%x_%j.log

path='/gpfs/commons/home/tlin/data'
clump_path='/gpfs/commons/home/tlin/output/cT/bellenguez/qc_on_base'

awk 'NR!=1{print $3}' $clump_path/bellenguez_clump_chr${chr}.clumped > $clump_path/chr${chr}.valid.snp

~/plink \
--bfile /gpfs/commons/home/tlin/data/biallelic/${chr}_filt \
--score /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/bellenguez_2021/bellenguez_2021_final_qc.tsv 3 4 7 header \
--q-score-range range_list.txt /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/bellenguez_2021/bellenguez_SNP.pvalue \
--extract $clump_path/chr${chr}.valid.snp \
--out /gpfs/commons/home/tlin/output/cT/bellenguez/qc_on_base/correct_column/bellenguez_pT_chr${chr}



##kunkle
#clump_path='/gpfs/commons/home/tlin/output/cT/kunkle/'
#clump_path='/gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap/clumping'
#--score ${path}/kunkle_etal_Stage1_info.tsv 3 4 8 header \
#--q-score-range range_list.txt $path/kunkle_check.pvalue \
#--out /gpfs/commons/home/tlin/output/cT/kunkle/kunkle_pT_chr${chr}


##$3 for SNP_ID,$4 for effective allele info, $7 for effect size estimate 

##QCed
#--score /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/bellenguez_2021/bellenguez_2021_final_qc.tsv 3 4 9 header \
#--q-score-range range_list.txt /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/bellenguez_2021/bellenguez_SNP.pvalue \



#no_qc
#--bfile /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/plink_biallelic/qc/ADSP_qc_chr${chr} \
#--score /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/bellenguez_2021/bellenguez_2021_final_rename.tsv 3 4 9 header \
#--q-score-range range_list.txt /gpfs/commons/home/tlin/data/bellenguez.pvalue \

