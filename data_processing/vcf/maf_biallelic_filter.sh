#!/bin/bash
#SBATCH --job-name=vcf_filt
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=70G
#SBATCH --time=35:00:00
#SBATCH --output=/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_preview/annotated_chunk_biallelic/%x_%j.log 

## the lower one doesnt work. I ended up using remove_non_biallelic.sh
source /gpfs/commons/groups/knowles_lab/software/anaconda3/bin/activate
conda activate polyfun

#path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/correct_chr_vcf'
#bcftools view -m2 -M2 -v snps -i 'MAF > 0.01' ADSP_annotated_chr${chr}.correct.vcf.gz -Ov -o filt/ADSP_annotated_chr${chr}.filt.vcf
#vcftools --vcf $path/filt/ADSP_annotated_chr${chr}.filt.vcf --max-missing 0.01 --min-alleles 2 --max-alleles 2 --recode --stdout > $path/filt/ADSP_annotated_chr${chr}.genofilt.vcf

path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_preview/'
bgzip -d $path/annotated_chunk/ADSP.chr${chr}.chunk${i}.vcf.bgz
vcftools --vcf $path/annotated_chunk/ADSP.chr${chr}.chunk${i}.vcf --max-missing 0.01 --min-alleles 2 --max-alleles 2 --recode --stdout > /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_preview/annotated_chunk_biallelic/ADSP.chr${chr}.chunk${i}.vcf