merge='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_preview/plink_hg38/qc/36k_ADSP_merge_list.txt'
 
echo /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_preview/plink_hg38/qc/ADSP.chr1.chunk1 > $merge

for chr in {1..5}
do
    chunk_num=$(ls /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_preview/plink_hg38| grep chr"${chr}".chunk | grep bed | wc -l)

    for ((i=1; i<=$chunk_num; i++)); 
        do
            echo /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_preview/plink_hg38/qc/ADSP.chr${chr}.chunk${i} >> $merge
        done
done



tail -n+3 "$merge" > "$merge.tmp" && mv "$merge.tmp" "$merge"
