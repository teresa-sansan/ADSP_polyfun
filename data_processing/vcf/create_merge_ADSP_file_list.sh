merge='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/plink/keep_biallelic/ADSP_merge_list.txt'
touch $merge
for i in {2..22}
do
#echo /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/plink/ADSP_annotated_chr${i} >>  /gpfs/commons/home/tlin/data/ADSP_merge_list.txt 
echo /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/plink/keep_biallelic/chr${i}.uni >> $merge
done
