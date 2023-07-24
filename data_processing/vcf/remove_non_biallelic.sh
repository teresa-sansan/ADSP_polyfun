#!/bin/bash
#SBATCH --job-name=bgzip
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=30G
#SBATCH --time=56:00:00
#SBATCH --output=/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_preview/annotated_chunk_biallelic/%x_%j.log

source /gpfs/commons/groups/knowles_lab/software/anaconda3/bin/activate
conda activate polyfun 
   
echo submit chunk num $i
bcftools view -m2 -M2 -v snps /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_preview/annotated_chunk/ADSP.chr${chr}.chunk${i}.vcf.bgz | bcftools filter -e 'F_MISSING > 0.1 || MAF <= 0.001' | bgzip -c >  /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_preview/annotated_chunk_biallelic/ADSP.chr${chr}.chunk${i}.vcf.bgzip




## a wrap script that works perfectly for the 36k bigwig

# for chr in {10..22}
# do
#     chunk_num=$(ls /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_preview/plink_hg38| grep chr"${chr}".chunk | grep bed | wc -l)
#     for ((i=1; i<=$chunk_num; i++)); 
#         do
#             echo check annotated_chunk_biallelic/ADSP.chr${chr}.chunk${i}
#             #sbatch --mem=50G -c 10 --wrap="bcftools view /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_preview/annotated_chunk/ADSP.chr${chr}.chunk${i}.vcf.bgz --regions chr${i}  > /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_preview/annotated_chunk_biallelic/ADSP.chr${chr}.chunk${i}"
#             #sbatch --mem=50G -c 10 --wrap="bcftools view -m2 -M2 -v snps /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_preview/annotated_chunk/ADSP.chr${chr}.chunk${i}.vcf.bgz | bgzip -c > /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_preview/annotated_chunk_biallelic/ADSP.chr${chr}.chunk${i}.vcf.gz"
#             sbatch --mem=50G -c 10 --wrap="source /gpfs/commons/groups/knowles_lab/software/anaconda3/bin/activate && conda activate polyfun && bgzip -c /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_preview/annotated_chunk_biallelic/ADSP.chr${chr}.chunk${i}  > /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_preview/annotated_chunk_biallelic/ADSP.chr${chr}.chunk${i}.vcf.bgz"
         
#         done
# done

# ##11-22