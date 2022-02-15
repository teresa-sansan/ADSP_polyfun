
for i in {15..22}
do
~/plink \
  --bfile /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/plink_biallelic/ADSP_chr$i \
  --maf 0.01 \
  --geno 0.01 \
  --mind 0.01 \
  --make-bed \
  --out /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/plink_biallelic/qc/ADSP_chr$i

done

#  --hwe 1e-6 \ ## Didn't include. Because we dont know which sample is control.
