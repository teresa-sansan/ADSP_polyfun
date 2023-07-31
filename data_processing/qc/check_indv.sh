
#cd /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_preview/plink_hg38_qc/processing_files
cd /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_preview/plink_hg38_qc/26k_processed
pwd
for chr in {20..22}
do  
    path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_preview/'
    chunk_num=$(ls $path/plink_hg38| grep chr"${chr}".chunk | grep bed | wc -l)
    echo chr $chr
    for ((i=1; i<=$chunk_num; i++)); 
        do 
        num=$(cat ADSP.mind.chr${chr}.chunk${i}.imiss|tail -n+2| tr -s ' ' | cut -d ' ' -f 6| uniq)
        echo chunk $i, $num
        done
    echo
done

## variant_count.sh

