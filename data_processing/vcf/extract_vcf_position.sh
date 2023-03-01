#!/bin/bash
#SBATCH --job-name=extract_position
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=15G
#SBATCH --time=22:00:00
#SBATCH --output=/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/only_position/%x%j.log

cd /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37

zcat ADSP_annotated_chr${chr}.vcf.gz| grep -v '##'| cut -f 1-10 > only_position/ADSP_annotated_hg37_chr${chr}.vcf.tsv


