#!/bin/bash
#SBATCH --job-name=check_correct_chr
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=150G
#SBATCH --time=22:00:00
#SBATCH --output=/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/correct_chr_vcf/%x%j.log

## this script is for extracting SNPs from the correct CHR. (because liftover tends to make it mismap)
cd /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37
for i in 4 6 
do
echo start chr$i
sbatch --mem=15G -c 10 --wrap="bcftools view ADSP_annotated_chr${i}.vcf.gz --regions chr${i} > /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/correct_chr_vcf/ADSP_annotated_chr${i}.correct.vcf"
done


