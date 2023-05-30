#!/bin/bash
for chr in {11..22}
do
    path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_preview/'
    chunk_num=$(ls $path/plink_hg38| grep chr"${chr}".chunk | grep bed | wc -l)
    echo submit chr ${chr} 

     for ((i=1; i<=$chunk_num; i++)); 
     do
     echo submit chunk num $i
     sbatch --export=chr=$chr,i=$i  extract_rsid_from_vcf.sh
     done
    echo
done