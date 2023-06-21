#!/bin/bash
#SBATCH --job-name=filt_indv
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=5G
#SBATCH --time=22:00:00
#SBATCH --output=/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_preview/plink_hg38_qc/processing_files/%x%j.log


PLINK_DIR='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_preview/plink_hg38_qc/processing_files'

~/plink \
  --bfile $PLINK_DIR/ADSP.mind.chr${chr}.chunk${chunk} \
  --remove $PLINK_DIR/ind_to_remove.tsv \
  --make-bed --out $PLINK_DIR/ADSP.finqc.chr${chr}.chunk${chunk}