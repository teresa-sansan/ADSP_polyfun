#!/bin/bash
#SBATCH --job-name=pT
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=50G
#SBATCH --time=2:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/cT/bellenguez/%x_%j.log

path='/gpfs/commons/home/tlin/data'



##qc on base
clump_path='/gpfs/commons/home/tlin/output/cT/bellenguez/qc_on_base'
#awk 'NR!=1{print $3}' $clump_path/bellenguez_clump_chr${chr}.clumped > $clump_path/chr${chr}.valid.snp

~/plink \
--bfile /gpfs/commons/home/tlin/data/biallelic/${chr}_filt \
--score /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/bellenguez_2021/bellenguez_2021_final_qc.tsv 3 4 7 header \
--q-score-range range_list.txt /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/bellenguez_2021/bellenguez_SNP.pvalue \
--extract $clump_path/chr${chr}.valid.snp \
--out /gpfs/commons/home/tlin/output/cT/bellenguez/qc_on_base/bellenguez_pT_chr${chr}




##qc_on_target
clump_path='/gpfs/commons/home/tlin/output/cT/bellenguez/qc_on_target'

~/plink \
--bfile /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/plink_biallelic/qc/ADSP_qc_chr${chr} \
--score /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/bellenguez_2021/bellenguez_2021_final_rename.tsv 3 4 7 header \
--q-score-range range_list.txt /gpfs/commons/home/tlin/data/bellenguez.pvalue \
--extract $clump_path/chr${chr}.valid.snp \
--out /gpfs/commons/home/tlin/output/cT/bellenguez/qc_on_target/bellenguez_pT_chr${chr}



##$3 for SNP_ID,$4 for effective allele info, $7 for effect size estimate 

