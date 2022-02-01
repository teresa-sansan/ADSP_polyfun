#!/bin/bash
#SBATCH --job-name=PRSice_prs
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=100G
#SBATCH --time=15:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/prs/PRSice/%x_%j.log 


Rscript /gpfs/commons/home/tlin/PRSice/PRSice.R \
   --prsice ~/PRSice/PRSice_linux \
   --base /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/bellenguez_2021/bellenguez_2021_final.tsv.gz \
   --target /gpfs/commons/home/tlin/data/biallelic/#_filt  \
   --binary-target T \
   --no-regress \
   --stat EFFECT \
   --A1 ALLELE1 \
   --A2 ALLELE2 \
   --snp RSID \
   --pvalue PVALUE \
   --beta \
   --out ~/output/prs/PRSice




#   --pheno /gpfs/commons/home/tlin/data/ADSP_pheno_merge_beta.tsv \
