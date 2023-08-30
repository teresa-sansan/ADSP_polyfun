#!/bin/bash
#SBATCH --job-name=36k_plink_ibd
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=30G
#SBATCH --time=5:00:00
#SBATCH --output=/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_QC/annotated_hg37_plink_ibd/%x_%j.log

source /gpfs/commons/groups/knowles_lab/software/anaconda3/bin/activate
conda activate polyfun 
## keep = SUBJID from /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/phenotype_file/release_36K/pheno_redo_his.tsv
## note the maf is 0.1%

chr=22
~/plink2/plink2 \
--const-fid \
--keep /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/phenotype_file/release_36K/keep_indv_ADSP_IBD.txt \
--make-bed \
--out /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_QC/plink_ibd/ADSP.chr${chr} \
--vcf /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_QC/annotated/ADSP.chr${chr}.vcf.gz


#--maf 0.001 \


#plink --vcf "$1" --out "$OUT$2" --const-fid --maf 0.000000000000001 --make-bed --keep "$EXTRACT"keep.txt
#--allow-extra-chr