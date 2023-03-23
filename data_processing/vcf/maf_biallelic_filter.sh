#!/bin/bash
#SBATCH --job-name=vcf_filt
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=70G
#SBATCH --time=5:00:00
#SBATCH --output=/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/correct_chr_vcf/filt/%x_%j.log 


source /gpfs/commons/groups/knowles_lab/software/anaconda3/bin/activate
conda activate polyfun

path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/correct_chr_vcf'
#bcftools view -m2 -M2 -v snps -i 'MAF > 0.01' ADSP_annotated_chr${chr}.correct.vcf.gz -Ov -o filt/ADSP_annotated_chr${chr}.filt.vcf
vcftools  --vcf $path/filt/ADSP_annotated_chr${chr}.filt.vcf --max-missing 0.01 --recode --out $path/filt/ADSP_annotated_chr${chr}.genofilt.vcf