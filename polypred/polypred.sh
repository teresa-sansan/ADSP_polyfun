#!/bin/bash
#SBATCH --job-name=Polypred
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=100G
#SBATCH --time=15:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap/polypred/%x_%j.log 


#/gpfs/commons/home/tlin/output/bellenguez/bellenguez_all_2/polypred/adjust_beta/%x_%j.log


cd /gpfs/commons/home/tlin/polyfun_omer_repo
source /gpfs/commons/groups/knowles_lab/software/anaconda3/bin/activate 
conda activate polyfun
path='/gpfs/commons/home/tlin/output/bellenguez/bellenguez_all_2'

prefix='PLINKupdate'

python polypred.py \
	--predict \
	--betas ${path}/finemap_snpvar_constrained/max_snp_${max_snp}/aggregrate.all.txt.gz \
	--output-prefix /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap/polypred/update_PLINK_${max_snp} \
	--plink-exe ~/plink \
	/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/plink_biallelic/ADSP_chr*.bed


##updateRSID


#	--betas /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap/max_snp_${max_snp}/aggregrate.all.txt.gz \
 

##SbayesR
#	--beta /gpfs/commons/home/tlin/output/sbayesR/mergebeta.betas \
#       --output-prefix /gpfs/commons/home/tlin/output/sbayesR/polypred/polypred.pred \


