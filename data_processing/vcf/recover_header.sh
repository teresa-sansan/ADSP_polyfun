#!/bin/bash
#SBATCH --job-name=vcf_header
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=100G
#SBATCH --time=50:00:00
#SBATCH --output=/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/correct_chr_vcf/filt/nodot/%x_%j.log

path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/correct_chr_vcf/filt'
cd $path
zcat ADSP_annotated_chr${chr}.genofilt.vcf.gz |grep '#' > nodot/ADSP_annotated_no_dot_chr${chr}.genofilt.vcf
zcat ADSP_annotated_no_dot_chr${chr}.genofilt.vcf.gz >>  nodot/ADSP_annotated_no_dot_chr${chr}.genofilt.vcf