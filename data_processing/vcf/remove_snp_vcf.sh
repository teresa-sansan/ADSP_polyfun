#!/bin/bash
#SBATCH --job-name=remove_Dot
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=80G
#SBATCH --time=10:00:00
#SBATCH --output=/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/correct_chr_vcf/filt/%x%j.log

path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/correct_chr_vcf/filt'
zcat $path/ADSP_annotated_chr${chr}.genofilt.vcf.gz |  grep -v '#'| awk '{if ($3 != ".") print $0}'|gzip > $path/ADSP_annotated_no_dot_chr${chr}.genofilt.vcf.gz
