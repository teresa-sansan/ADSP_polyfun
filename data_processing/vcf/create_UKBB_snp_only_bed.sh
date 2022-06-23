#!/bin/sh 
#SBATCH --job-name=extract_UKBB_only
#SBATCH --mail-type=FAIL,END 
#SBATCH --mail-user=tlin@nygenome.org 
#SBATCH --mem=50G 
#SBATCH --time=5:00:00
#SBATCH --output=/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/plink/%x.log 


cd /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/plink

~/plink \
	--bfile ADSP_annotated_fixed \
	--extract /gpfs/commons/groups/knowles_lab/data/ldsc/polyfun/baselineLF2.2.UKB/rsid_only.txt \
	--make-bed --out ADSP_UKBB_only \

