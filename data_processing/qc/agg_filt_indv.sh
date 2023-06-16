ind_path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_preview/plink_hg38_qc'
output_name='ind_removed.tsv'
echo 'ID CHR CHUNK' > $ind_path/$output_name

for chr in {1..22}
do
    path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_preview/'
    chunk_num=$(ls $path/plink_hg38| grep chr"${chr}".chunk | grep bed | wc -l)

     for ((i=1; i<=$chunk_num; i++)); 
     do
        if [ -f $ind_path/ADSP.mind.chr${chr}.chunk${i}.irem ]; then
            cat $ind_path/ADSP.mind.chr${chr}.chunk${i}.irem| cut -f 2 | awk -v chr="$chr" -v chunk="$i" '{ print $0, chr, chunk }'>> $ind_path/$output_name
        fi  
     done

done    
echo Done aggregating, now counting the occurances...

awk '{count[$1]++} END{for (val in count) print val, count[val]}' $ind_path/$output_name > $ind_path/ind_count.txt
echo Done counting

cat ind_count.txt| awk -F ' ' '{ if ($2 > 6) print $1}' > ind_to_remove.txt
echo extracting the IDs to remove