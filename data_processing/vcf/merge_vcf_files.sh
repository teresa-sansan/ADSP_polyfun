#!/bin/bash
#SBATCH --job-name=merge_vcf_files
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=900G 
#SBATCH --time=15:00:00
#SBATCH --partition bigmem
#SBATCH --output=/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/correct_chr_vcf/filt/%x%j.log


cd /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/correct_chr_vcf/filt/
bcftools concat  -a -f concat_list.txt  -o ADSP_annotated_merged.vcf
