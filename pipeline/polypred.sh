#!/bin/bash
#SBATCH --job-name=Polypred
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=100G
#SBATCH --time=15:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/bellenguez/bellenguez_all_2/polypred/%x_%j.log


cd /gpfs/commons/home/tlin/polyfun_omer_repo
source /gpfs/commons/groups/knowles_lab/software/anaconda3/bin/activate 
conda activate polyfun
python polypred.py \
	--predict \
	--betas /gpfs/commons/home/tlin/output/bellenguez/bellenguez_all_2/finemap_snpvar_constrained/max_snp_${max_snp}/aggregrate.all.txt.gz \
	--output-prefix /gpfs/commons/home/tlin/output/bellenguez/bellenguez_all_2/polypred/polypred_bellenguez_max${max_snp} \
	--plink-exe ~/plink \
	/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/plink_biallelic/ADSP_chr*.bed

#	--betas /gpfs/commons/home/tlin/output/kunkle_all_2/finemap/max_snp_1/chr${chr}.aggregrate.all.txt.gz \
#	/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/filtered_data/ADSP_chr*.bed
# 	--output-prefix /gpfs/commons/home/tlin/output/kunkle_all/polypred/polypred_kunkle_all2_chr${chr} \

