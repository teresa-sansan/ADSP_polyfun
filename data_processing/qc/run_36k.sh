PLINK_DIR='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_preview/plink_hg38_qc'
#touch >> $PLINK_DIR/qc_indv.count.txt
for chr in {1..19}
do
    path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_preview/'
    chunk_num=$(ls $path/plink_hg38| grep chr"${chr}".chunk | grep bed | wc -l)
    echo chr ${chr} 

    for ((i=1; i<=$chunk_num; i++)); 
    do
      #cat $PLINK_DIR/ADSP.finqc.chr${chr}.chunk${i}.log| grep "people pass filters and QC" >> $PLINK_DIR/qc_indv.count.txt
      sbatch --export=chr=$chr,chunk=$i filt_maf_plink.sh
    done
done