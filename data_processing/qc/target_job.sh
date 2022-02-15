#!/bin/bash
#SBATCH --job-name=target_QC
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=150G
#SBATCH --time=3:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/others/qc/%x_%j.log

## Filtered on MAF, genotype missingness rate

~/plink \
  --bfile /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/plink_biallelic/ADSP_chr$i \
  --maf 0.01 \
  --geno 0.01 \
  --mind 0.01 \
  --make-bed \
  --out /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/plink_biallelic/qc/ADSP_qc_chr$i

#  --hwe 1e-6 \ ## Didn't include. Because we dont know which sample is control.


## Remove highly correlated SNPs
~/plink \
  --bfile /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/plink_biallelic/qc/ADSP_qc_chr$i \
  --indep-pairwise 200 50 0.25 \
  --out /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/plink_biallelic/qc/ADSP_qc_chr$i


## Calculate hetrozygosity rates
~/plink \
  --bfile /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/plink_biallelic/qc/ADSP_qc_chr$i \
  --extract /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/plink_biallelic/qc/ADSP_qc_chr$i.prune.in \
  --keep /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/plink_biallelic/qc/ADSP_qc_chr$i.fam \
  --het \
  --out /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/plink_biallelic/qc/ADSP_qc_chr$i

