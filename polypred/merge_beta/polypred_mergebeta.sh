#!/bin/bash
#SBATCH --job-name=Polypred_mergebeta
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=100G
#SBATCH --time=15:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/sbayesR/polypred/%x_%j.log

cd /gpfs/commons/home/tlin/polyfun_omer_repo/re-clone
source /gpfs/commons/groups/knowles_lab/software/anaconda3/bin/activate
conda activate polyfun


## only use one chr to see the result (faster debug)
python polypred.py --combine-betas \
  --betas /gpfs/commons/home/tlin/output/bellenguez/bellenguez_all_2/finemap_snpvar_constrained/max_snp_10/aggregrate.all.txt.gz,~/output/sbayesR/bellenguez_agg_snpRes \
  --pheno /gpfs/commons/home/tlin/data/ADSP_pheno_merge_beta.tsv \
  --output-prefix /gpfs/commons/home/tlin/output/sbayesR/polypred/updated_polypred/0204_pull \
  --plink-exe ~/plink /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/plink_biallelic/ADSP_chr*.bed

#  --betas /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap/max_snp_10/aggregrate.all.txt.gz,~/output/sbayesR/bellenguez_agg_snpRes \
#  --betas ~/output/sbayesR/bellenguez_agg_snpRes,/gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap/max_snp_10/aggregrate.all.txt.gz \
