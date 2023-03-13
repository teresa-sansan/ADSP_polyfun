#!/bin/bash
#SBATCH --job-name=gzip
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=150G
#SBATCH --time=22:00:00
#SBATCH --output=/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/correct_chr_vcf/%x%j.log


cd /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/correct_chr_vcf
for i in {7..22}
do
echo start chr$i
sbatch --mem=15G -c 10 --wrap="gzip ADSP_annotated_chr${i}.correct.vcf"
done


