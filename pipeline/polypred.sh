#!/bin/bash
#SBATCH --job-name=Polypred_susie
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=100G
#SBATCH --time=15:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/prs/bellenguez_all_2/clump/%x_%j.log 


#/gpfs/commons/home/tlin/output/bellenguez/bellenguez_all_2/polypred/adjust_beta/%x_%j.log


cd /gpfs/commons/home/tlin/polyfun_omer_repo
source /gpfs/commons/groups/knowles_lab/software/anaconda3/bin/activate 
conda activate polyfun
path='/gpfs/commons/home/tlin/output/bellenguez/bellenguez_all_2'

prefix='polypred_adj_beta_bellenguez_all2'

python polypred_PT.py \
	--predict \
	--betas ${path}/finemap_snpvar_constrained/max_snp_${max_snp}/aggregrate_adjust_beta.all.txt.gz \
	--output-prefix /gpfs/commons/home/tlin/output/prs/bellenguez_all_2/clump/test_PT \
	--plink-exe ~/plink \
	/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/plink_biallelic/ADSP_chr*.bed 
#	--output-prefix /gpfs/commons/home/tlin/output/bellenguez/bellenguez_all_2/polypred/adjust_beta \
#	--output-prefix ${path}/polypred/adjust_beta/${prefix}_${max_snp} \

#path='/gpfs/commons/home/tlin/output/bellenguez/bellenguez_bl'
#prefix='polypred_bellenguez_susie'

#path='/gpfs/commons/home/tlin/output/bellenguez/bellenguez_all_2/finemap_snpvar_constrained/'

#	--betas /gpfs/commons/home/tlin/output/kunkle_all_2/finemap/max_snp_1/chr${chr}.aggregrate.all.txt.gz \
#	/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/filtered_data/ADSP_chr*.bed
# 	--output-prefix /gpfs/commons/home/tlin/output/kunkle_all/polypred/polypred_kunkle_all2_chr${chr} \

