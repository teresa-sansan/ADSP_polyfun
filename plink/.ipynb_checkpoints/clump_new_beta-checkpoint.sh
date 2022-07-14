#!/bin/bash
#SBATCH --job-name=clump_new_beta
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=40G
#SBATCH --time=10:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/cT/chr_sep_plink/wightman/new_beta/%x_%j.log

## Even though this step doesn't need beta, but just to make sure they are the same. 

sumstats='/gpfs/commons/home/tlin/output/wightman/fixed_0224/finemap/max_snp_10/agg_all_new_beta.tsv' 

## Pvalue, SNP

##qc 
#--bfile /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/plink_biallelic/qc/ADSP_qc_chr${chr} \

~/plink \
--bfile /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/plink_biallelic/filt/${chr}_filt \
--clump-p1 1 \
--clump-r2 0.1  \
--clump-kb 250  \
--clump $sumstats \
--clump-snp-field SNP \
--clump-field P \
--out /gpfs/commons/home/tlin/output/cT/chr_sep_plink/wightman/new_beta/wightman_${chr}



## kunkle
#sumstats='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/new_beta/kunkle_max_snp_10.aggregate.tsv.gz'
##/gpfs/commons/home/tlin/output/cT/kunkle/fixed_0224/new_beta/kunkle_${chr}

