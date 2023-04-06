#!/bin/bash
#SBATCH --job-name=bgzip
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=100G
#SBATCH --time=40:00:00
#SBATCH --output=/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/plink/vcf_filt/%x_%j.log



path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/correct_chr_vcf/filt'
echo start chr${i}
source /gpfs/commons/groups/knowles_lab/software/anaconda3/bin/activate
conda activate polyfun
gunzip -c ${path}/ADSP_annotated_chr${i}.genofilt.vcf.gz | bgzip > ${path}/bgzip/ADSP_annotated_chr${i}.genofilt.vcf.bgz
#sbatch --mem=90G -c 5 --wrap="gzip ${path}/ADSP_annotated_chr${i}.genofilt.vcf > ${path}/gzip/ADSP_annotated_chr${i}.genofilt.vcf.gz"
