#!/bin/bash
#SBATCH --job-name=bgzip
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=15G
#SBATCH --time=16:00:00
#SBATCH --output=/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_preview/annotated_chunk_biallelic/%x_%j.log

source /gpfs/commons/groups/knowles_lab/software/anaconda3/bin/activate
conda activate polyfun 
for chr in {20..22}
do
    path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_preview/'
    chunk_num=$(ls $path/plink_hg38| grep chr"${chr}".chunk | grep bed | wc -l)
    echo submit chr ${chr} 

     for ((i=1; i<=$chunk_num; i++)); 
     do
        sbatch --export=chr=$chr,chunk=$i maf_biallelic_filter.sh
    #  echo submit chunk num $i
    #  bgzip -c /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_preview/annotated_chunk_biallelic/ADSP.chr${chr}.chunk${i}  > /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_preview/annotated_chunk_biallelic/ADSP.chr${chr}.chunk${i}.vcf.bgzip
    #  rm /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_preview/annotated_chunk_biallelic/ADSP.chr${chr}.chunk${i}
         
     #sbatch --export=chr=$chr,i=$i  extract_rsid_from_vcf.sh
     #sbatch --export=chr=$chr,chunk=$i filt_maf_plink.sh
     #sbatch --export=chr=$chr,chunk=$i maf_biallelic_filter.sh
     done
    echo
done