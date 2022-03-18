touch /gpfs/commons/home/tlin/data/ADSP_qc_merge_list.txt
for i in {2..22}
do
#echo /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/plink_biallelic/filt/${i}_filt >>  /gpfs/commons/home/tlin/data/ADSP_merge_list.txt 
echo /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/plink_biallelic/qc/ADSP_qc_chr${i} >>  /gpfs/commons/home/tlin/data/ADSP_qc_merge_list.txt 
done
