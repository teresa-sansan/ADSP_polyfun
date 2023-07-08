merge='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_preview/plink_hg38_qc/36k_ADSP_merge_list.txt' 
echo /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_preview/plink_hg38_qc/processing_files/rm_multi_alleleic/ADSP.chr1.chunk1 > $merge

for chr in {1..22}
do
    chunk_num=$(ls /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_preview/plink_hg38| grep chr"${chr}".chunk | grep bed | wc -l)

    for ((i=1; i<=$chunk_num; i++)); 
        do
            echo /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_preview/plink_hg38_qc/processing_files/rm_multi_alleleic/ADSP.chr${chr}.chunk${i} >> $merge
        done
done


tail -n+2 "$merge" > "$merge.tmp" && mv "$merge.tmp" "$merge"
