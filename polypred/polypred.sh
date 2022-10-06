#!/bin/bash
#SBATCH --job-name=Polypred_new_wightman
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=50G
#SBATCH --time=5:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/wightman/wightman_check_1003/all_anno/finemap/polypred/%x_%j.log 

cd /gpfs/commons/home/tlin/polyfun_omer_repo
source /gpfs/commons/groups/knowles_lab/software/anaconda3/bin/activate 
conda activate polyfun

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

#path='/gpfs/commons/home/tlin/output/kunkle/kunkle_fixed_0224_annotations/new_susie'

## bellenguez
## susie
#path='/gpfs/commons/home/tlin/output/bellenguez/bellenguez_fixed_0224_annotations'
#path='/gpfs/commons/home/tlin/output/bellenguez/old/bellenguez_fixed_0224'

##wightman

#path=/gpfs/commons/home/tlin/output/wightman/fixed_0224/finemap
#	--output-prefix $path/polypred/max_snp_${max_snp}.new_beta_wightman_polypred.tsv \
#	--betas $path/max_snp_${max_snp}/agg_all_new_beta.tsv.gz \
#	--output-prefix $path/polypred/max_snp_${max_snp}_new_beta_polypred.tsv \
#path='/gpfs/commons/home/tlin/output/wightman/fixed_0224/susie/finemap_fixed_assertion_susie_iter'
#path='/gpfs/commons/home/tlin/output/wightman/fixed_0224_annotations'

#path='/gpfs/commons/home/tlin/output/wightman/wightman_fixed_0224/finemap'
#python ~/polyfun_omer_repo/polypred_new_beta.py \
path='/gpfs/commons/home/tlin/output/wightman/wightman_check_1003/all_anno/finemap/'


## jansen
path='/gpfs/commons/home/tlin/output/jansen/finemap'
#path='/gpfs/commons/home/tlin/output/jansen/susie'

python polypred.py \
	--predict \
	--betas $path/max_snp_${max_snp}/aggregate.all.txt \
	--output-prefix $path/polypred/max_snp_${max_snp}_polypred.tsv \
	--plink-exe ~/plink \
	/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/plink/ADSP_qc_all/ADSP_qc_all_*.bed


#bellenguez
#	--betas $path/finemap/max_snp_${max_snp}/aggregate.all.txt \

#$path/${anno}/max_snp_${max_snp}/aggregate.all.txt \
#	--betas /gpfs/commons/home/tlin/output/wightman/wightman_fixed_0224/finemap/max_snp_3/aggregate.all.txt.gz \
#	--output-prefix $path/polypred_new_plink/max_snp_${max_snp}_polypred.tsv \

##aggregate.all.txt.gz 
