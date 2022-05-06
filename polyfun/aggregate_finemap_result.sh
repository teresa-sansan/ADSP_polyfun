#!/bin/bash
#SBATCH --job-name=aggregrate_test_rescue_IBSS
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=100G
#SBATCH --time=3:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/bellenguez/bellenguez_fixed_0224/finemap/try_rescue_not_converge/%x_%j.log

cd /gpfs/commons/home/tlin/polyfun_omer_repo

source /gpfs/commons/groups/knowles_lab/software/anaconda3/bin/activate
conda activate polyfun


##bellenguez
#bellenguez='/gpfs/commons/home/tlin/data/bellenguez_munged.parquet'
bellenguez='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/Bellenguez_et_al_2021_hg37.munged.parquet'


##bellenguez_updateRSID
bellenguez='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/bellenguez_2021/bellenguez_munged.parquet'

prefix='finemap_max_snp_3'
path='/gpfs/commons/home/tlin/output/bellenguez/bellenguez_fixed_0224/finemap/try_rescue_not_converge'
#prefix='finemap_bellenguez'
#path='/gpfs/commons/home/tlin/output/bellenguez/bellenguez_fixed_0224/finemap'
#path='/gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap'

##wightman
#wightman='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/Wightman_et_al_2021_hg37_ldsc.tsv.gz'
#wightman='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/wightman_2021/wightman_2021_fixed.parquet'
#prefix='finemap_wightman'
#path='/gpfs/commons/home/tlin/output/wightman/finemap'

#kunkle='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/kunkle_2019/AD_Kunkle_etal_Stage1.parquet'
kunkle='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/kunkle_2019/Kunkle_etal_Stage1_qc.gz'

if false; then
for chr in {1..22}
do
	#python aggregate_finemapper_results_modified.py \
	python aggregate_finemapper_results_min.py \
		--out-prefix $path/max_snp_${max_snp}/$prefix \
		--sumstats $bellenguez \
		--out $path/max_snp_${max_snp}/min_chr${chr}.aggregrate.all.txt.gz \
		--allow-missing \
       		--chr $chr

done

fi



## this part is trying to rescue the region that fails to converge. 

if true; then

python aggregate_finemapper_results_fixed_convergence.py \
                --out-prefix $path/$prefix \
                --sumstats $bellenguez \
                --allow-missing \
                --out $path/aggregrate.all.txt.gz

fi
