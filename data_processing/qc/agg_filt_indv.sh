#ind_path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_preview/plink_hg38_qc/processing_files'
##26k
ind_path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_preview/plink_hg38_qc/26k_processed'
output_name='ind_removed.tsv'
echo 'ID CHR CHUNK' > $ind_path/$output_name

for chr in {1..22}
do
    path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_preview/'
    chunk_num=$(ls $path/plink_hg38| grep chr"${chr}".chunk | grep bed | wc -l)

     for ((i=1; i<=$chunk_num; i++)); 
     do
        if [ -f $ind_path/ADSP.mind.chr${chr}.chunk${i}.irem ]; then
            #cat $ind_path/ADSP.mind.chr${chr}.chunk${i}.irem| cut -f 2 | awk -v chr="$chr" -v chunk="$i" '{ print $0, chr, chunk }'>> $ind_path/$output_name
            cat $ind_path/ADSP.mind.chr${chr}.chunk${i}.imiss|tr -s ' '|tr ' ' '\t'| cut -f 3,5,6  > $ind_path/var_count_chr${chr}.chunk${i}.tsv 
        fi  
     done
done    
echo "Done calculating the missingness per indv per chunk, now summing up the number of total missing variants per individual..."
#awk '{ sum2[$1]+=$2} END { for (user in sum2) print user, sum2[user] }' $ind_path/var_count_chr* > $ind_path/check.tsv
awk '{sum2[$1]+=$2; sum3[$1]+=$3 } END { for (user in sum2) print user, sum2[user], sum3[user] }' $ind_path/var_count_chr* > $ind_path/check.tsv


#$ind_path/ind_missingvar_count.tsv
echo Done counting, saved at  $ind_path/check.tsv


# cat ind_count.txt| awk -F ' ' '{ if ($2 > 6) print $1}' > ind_to_remove.txt
# echo extracting the IDs to remove