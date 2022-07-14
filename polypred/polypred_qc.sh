#!/bin/bash
#SBATCH --job-name=Polypred_bellenguez_plinkUKBB
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=100G
#SBATCH --time=15:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/bellenguez/bellenguez_fixed_0224/finemap/try_rescue_not_converge/polypred/%x_%j.log 

cd /gpfs/commons/home/tlin/polyfun_omer_repo
source /gpfs/commons/groups/knowles_lab/software/anaconda3/bin/activate 
conda activate polyfun
path='/gpfs/commons/home/tlin/output/bellenguez/bellenguez_fixed_0224'

max_snp=10
python polypred.py \
	--predict \
	--betas /gpfs/commons/home/tlin/output/bellenguez/bellenguez_fixed_0224/finemap/try_rescue_not_converge/finemap_genomewide.tsv \
	--output-prefix /gpfs/commons/home/tlin/output/bellenguez/bellenguez_fixed_0224/finemap/try_rescue_not_converge/polypred/fixed_0224_UKBBonly \
	--plink-exe ~/plink \
	/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/plink/ADSP_UKBB/ADSP_UKBB_only_*.bed

	
##plink
#	/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/plink_biallelic/qc_maf_0.1/ADSP_qc_chr*.bed
#	/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/plink/ADSP_UKBB/*_filt.bed


##no_qced target
#	/gpfs/commons/home/tlin/data/biallelic/ADSP_chr*.bed


##fixed 0224
#	--betas $path/finemap/max_snp_${max_snp}/aggregrate.all.txt.gz \
#	--output-prefix $path/finemap/polypred/max_snp_${max_snp} \


##updateRSID
#	--betas /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap/max_snp_${max_snp}/aggregrate.all.txt.gz \
#	--output-prefix /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap/polypred_qc/max_snp_${max_snp} \
 

##SbayesR
#	--beta /gpfs/commons/home/tlin/output/sbayesR/mergebeta.betas \
#       --output-prefix /gpfs/commons/home/tlin/output/sbayesR/polypred/polypred.pred \


##updatePlink
#	--betas ${path}/finemap_snpvar_constrained/max_snp_${max_snp}/aggregrate.all.txt.gz \

