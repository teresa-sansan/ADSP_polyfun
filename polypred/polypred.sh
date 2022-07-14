#!/bin/bash
#SBATCH --job-name=Polypred
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=50G
#SBATCH --time=15:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/wightman/fixed_0224/susie/finemap_fixed_assertion_susie_iter/max_snp_10/%x_%j.log 

cd /gpfs/commons/home/tlin/polyfun_omer_repo
source /gpfs/commons/groups/knowles_lab/software/anaconda3/bin/activate 
conda activate polyfun

#max_snp=10
#path='/gpfs/commons/home/tlin/output/bellenguez/bellenguez_fixed_0224/finemap/try_rescue_not_converge/'
#path='/gpfs/commons/home/tlin/output/bellenguez/bellenguez_fixed_0224'
#path='/gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID'
#	--output-prefix /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap/polypred_qc/${max_snp} \
#	--betas $path/max_snp_${max_snp}/aggregrate.all.txt.gz \


##kunkle
##susie
#path='/gpfs/commons/home/tlin/output/kunkle/kunkle_fixed_0224/susie_finemap'
#path='/gpfs/commons/home/tlin/output/kunkle/kunkle_fixed_0224/finemap'
#output='	--output-prefix /gpfs/commons/home/tlin/output/bellen'
#	--betas $path/max_snp_${max_snp}/agg_kunkle_extract_1e-3.tsv \



##wightman
#path=/gpfs/commons/home/tlin/output/wightman/fixed_0224/finemap
#	--output-prefix $path/polypred/max_snp_${max_snp}.new_beta_wightman_polypred.tsv \
#	--betas $path/max_snp_${max_snp}/agg_all_new_beta.tsv.gz \
#	--output-prefix $path/polypred/max_snp_${max_snp}_new_beta_polypred.tsv \

path='/gpfs/commons/home/tlin/output/wightman/fixed_0224/susie/finemap_fixed_assertion_susie_iter'


max_snp=10

#python ~/polyfun_omer_repo/polypred_new_beta.py \

python polypred.py \
	--predict \
	--betas $path/max_snp_${max_snp}/aggregate.all.txt.gz \
	--output-prefix $path/polypred/max_snp_${max_snp}_polypred.tsv \
	--plink-exe ~/plink \
	/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/plink_biallelic/filt/*_filt.bed


##updateRSID
#	--betas /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap/max_snp_${max_snp}/aggregrate.all.txt.gz \
 

##SbayesR
#	--beta /gpfs/commons/home/tlin/output/sbayesR/mergebeta.betas \
#       --output-prefix /gpfs/commons/home/tlin/output/sbayesR/polypred/polypred.pred \


##updatePlink
#	--betas ${path}/finemap_snpvar_constrained/max_snp_${max_snp}/aggregrate.all.txt.gz \

