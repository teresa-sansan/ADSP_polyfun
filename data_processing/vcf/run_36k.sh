#!/bin/bash
#SBATCH --job-name=36k_plink
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=50G
#SBATCH --time=100:00:00
#SBATCH --output=/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_preview/plink_hg38/all/%x_%j.log

source /gpfs/commons/groups/knowles_lab/software/anaconda3/bin/activate
conda activate polyfun
for chr in {1..21}
do
    # path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_preview/'
    # chunk_num=$(ls $path/plink_hg38| grep chr"${chr}".chunk | grep bed | wc -l)
    # echo submit chr ${chr} 

     #for ((i=1; i<=$chunk_num; i++)); 
     #do
        echo submit chunk num $chr
        #sbatch --export=chr=$chr,i=$i remove_non_biallelic.sh
        #sbatch --export=chr=$chr,chunk=$i create_plink.sh
       
    #  bgzip -c /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_preview/annotated_chunk_biallelic/ADSP.chr${chr}.chunk${i}  > /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_preview/annotated_chunk_biallelic/ADSP.chr${chr}.chunk${i}.vcf.bgzip
    #  rm /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_preview/annotated_chunk_biallelic/ADSP.chr${chr}.chunk${i}
         
     #sbatch --export=chr=$chr,i=$i  extract_rsid_from_vcf.sh
     #sbatch --export=chr=$chr,chunk=$i ../qc/target_qc_indv.sh
     #sbatch --export=chr=$chr,chunk=$i maf_biallelic_filter.sh
     #done

    #sbatch --export=i=$chr check_mismatch_liftover.sh
    sbatch --export=i=$chr extract_correct_vcf.sh
    echo
done