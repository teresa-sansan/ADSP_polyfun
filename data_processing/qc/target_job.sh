#!/bin/bash
#SBATCH --job-name=target_QC
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=150G
#SBATCH --time=3:00:00
#SBATCH --output=/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/plink/%x_%j.log

## Filtered on MAF, genotype missingness rate
path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/plink'

~/plink \
  --bfile $path/ADSP_annotated_fixed \
  --maf 0.001 \
  --geno 0.01 \
  --make-bed \
  --out $path/ADSP_qc_variant

#  --hwe 1e-6 \ ## Didn't include. Because we dont know which sample is control.
#  --mind 0.01 \ ## removed, because we want to test out how the performance is when only considerating variants. 



~/plink \
  --bfile $path/ADSP_annotated_fixed \
  --maf 0.001 \
  --geno 0.01 \
  --mind 0.01 \ 
  --make-bed \
  --out $path/ADSP_qc_target

## qc on individual alone
if false; then 
~/plink \
  --bfile /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/plink_biallelic/filt/${i}_filt \
  --mind 0.01 \
  --make-bed \
  --out /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/plink_biallelic/qc_on_individual/ADSP_qc_chr$i
fi


## Remove highly correlated SNPs
if false;then
~/plink \
  --bfile /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/plink_biallelic/qc_maf_0.1/ADSP_qc_chr$i \
  --indep-pairwise 200 50 0.25 \
  --out /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/plink_biallelic/qc_maf_0.1/ADSP_qc_chr$i
fi

## Calculate hetrozygosity rates
if false; then
~/plink \
  --bfile /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/plink_biallelic/qc_maf_0.1/ADSP_qc_chr$i \
  --extract /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/plink_biallelic/qc_maf_0.1/ADSP_qc_chr$i.prune.in \
  --keep /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/plink_biallelic/qc_maf_0.1/ADSP_qc_chr$i.fam \
  --het \
  --out /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/plink_biallelic/qc_maf_0.1/ADSP_qc_chr$i
fi
