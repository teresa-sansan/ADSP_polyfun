#!/bin/sh
#SBATCH --job-name=plink_vcf_chr6
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=50G
#SBATCH --time=10:00:00
#SBATCH --output=/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/plink/keep_biallelic/%x.log

cd /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/plink
/gpfs/commons/home/tlin/plink --bfile ADSP_annotated_chr6 --snps-only --exclude ADSP_annotated-merge.missnp --make-bed --out keep_biallelic/chr6.uni
