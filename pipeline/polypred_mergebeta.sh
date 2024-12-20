#!/bin/bash
#SBATCH --job-name=Polypred_mergebeta
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=100G
#SBATCH --time=15:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/sbayesR/polypred/%x_%j.log



cd /gpfs/commons/home/tlin/polyfun_omer_repo
source /gpfs/commons/groups/knowles_lab/software/anaconda3/bin/activate
conda activate polyfun



python polypred_modified.py \
  --combine-betas \
  --betas /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap/max_snp_10/aggregrate.all.txt.gz,~/output/sbayesR/bellenguez_agg_snpRes \
  --pheno /gpfs/commons/home/tlin/data/ADSP_pheno_merge_beta.tsv \
  --output-prefix /gpfs/commons/home/tlin/output/sbayesR/mergebeta \
  --plink-exe ~/plink /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/plink_biallelic/ADSP_chr*.bed

#  --betas ~/output/sbayesR/bellenguez_agg_snpRes,/gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap/max_snp_10/aggregrate.all.txt.gz \
