for chr in {1..22}
do
echo run chr $chr

~/plink \
--bfile /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/plink_biallelic/filt/${chr}_filt \
--score /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/bellenguez_2021/Bellenguez_qc_interested_SNP_maxPIP.tsv 2 4 12 header \
--out /gpfs/commons/home/tlin/output/prs/bellenguez/updateRSID/interested_SNP/maxPIP/updateRSID_interested_SNP_maxPIP_chr${chr}

~/plink \
--bfile /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/plink_biallelic/filt/${chr}_filt \
--score /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/bellenguez_2021/Bellenguez_qc_interested_SNP_minPIP.tsv 2 4 12 header \
--out /gpfs/commons/home/tlin/output/prs/bellenguez/updateRSID/interested_SNP/minPIP/updateRSID_interested_SNP_minPIP_chr${chr}

done


## SNP, effect allele, effect size

