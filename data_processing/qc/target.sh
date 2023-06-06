#!/bin/bash
#SBATCH --job-name=target_qc
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=150G
#SBATCH --time=22:00:00
#SBATCH --output=/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/plink/vcf_filt/%x_%j.log


~/plink \
  --bfile /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/plink/vcf_filt/ADSP_annotated_merged \
  --geno 0.01 \
  --make-bed --out /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/plink/vcf_filt/ADSP_annotated_merged_qc

#done

#  --hwe 1e-6 \ ## Didn't include. Because we dont know which sample is control.
#--mind 0.01 \  ## done in vcf
#--maf 0.001 \ ## done in vcftools (set to 0.1%)