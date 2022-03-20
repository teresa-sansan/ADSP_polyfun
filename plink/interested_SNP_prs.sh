for chr in 2 3 4 6 7 8 11 14 15 16 17 19
do
~/plink \
--bfile /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/plink_biallelic/filt/${chr}_filt \
--score /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/bellenguez_2021/Bellenguez_interested_SNP.tsv 2 4 12 header \
--out /gpfs/commons/home/tlin/output/prs/bellenguez/updateRSID_interested_SNP_chr${chr}

~/plink \
--bfile /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/plink_biallelic/qc/ADSP_qc_chr${chr} \
--score /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/bellenguez_2021/Bellenguez_qc_interested_SNP.tsv 2 4 12 header \
--out /gpfs/commons/home/tlin/output/prs/bellenguez/updateRSID_qc_interested_SNP_chr${chr}



done


## SNP, effect allele, effect size

