#!/bin/bash
#SBATCH --job-name=target_qc_step1
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=150G
#SBATCH --time=22:00:00
#SBATCH --output=/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/plink_biallelic/check/%x_%j.log

for i in {1..22}
do
~/plink \
  --bfile /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/plink_biallelic/ADSP_chr$i \
  --maf 0.01 \
  --geno 0.01 \
  --mind 0.01 \
  --make-bed \
  --out /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/plink_biallelic/check/ADSP_chr$i

done

#  --hwe 1e-6 \ ## Didn't include. Because we dont know which sample is control.