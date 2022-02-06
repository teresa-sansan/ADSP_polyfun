#!/bin/bash
#SBATCH --job-name=PRSice_prs_wpheno_wcov
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=200G
#SBATCH --time=15:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/prs/PRSice/%x_%j.log 


Rscript /gpfs/commons/home/tlin/PRSice/PRSice.R \
   --prsice ~/PRSice/PRSice_linux \
   --base /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/bellenguez_2021/bellenguez_2021_final_rename.tsv.gz \
   --target /gpfs/commons/home/tlin/data/biallelic/#_filt  \
   --binary-target T \
   --pheno /gpfs/commons/home/tlin/data/ADSP_pheno_merge_beta.tsv \
   --ignore-fid \
   --pheno-col AD_status_final \
   --cov /gpfs/commons/home/tlin/script/polypred/merge_beta/all_phenotypes_cov.tsv \
   --stat EFFECT \
   --A1 A1 \
   --A2 A2 \
   --snp SNP \
   --pvalue PVALUE \
   --beta \
   --out ~/output/prs/PRSice/PRSice_pheno_cov



#   --pheno /gpfs/commons/home/tlin/data/ADSP_pheno_merge_beta.tsv \



