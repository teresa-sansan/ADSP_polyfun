#!/bin/bash
#SBATCH --job-name=aggregrate_result
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=100G
#SBATCH --time=3:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap/%x_%j.log



cd /gpfs/commons/home/tlin/polyfun_omer_repo

source /gpfs/commons/groups/knowles_lab/software/anaconda3/bin/activate
conda activate polyfun

kunkle='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/kunkle_2019/AD_Kunkle_etal_Stage1.parquet'
#bellenguez='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/bellenguez_2021/Bellenguez_2021_stage1.parquet'
bellenguez='/gpfs/commons/home/tlin/data/bellenguez_2021_final.tsv.gz'
path='/gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap'
prefix='finemap_bellenguez'


for chr in {1..22}
do
	python aggregate_finemapper_results_modified.py \
		--out-prefix $path/max_snp_${max_snp}/$prefix \
		--sumstats $bellenguez \
		--out $path/max_snp_${max_snp}/chr${chr}.aggregrate.all.txt.gz \
		--allow-missing \
       		--chr $chr 
done


#prefix='/gpfs/commons/home/tlin/output/kunkle_all/finemap/all_anno'
#output='~/output/kunkle_all/finemap/aggregrate.all.txt.gz'
#prefix='/gpfs/commons/home/tlin/output/bellenguez/bellenguez_all_2/finemap/max_snp_1/finemap_bellenguez_all_2'
#path='/gpfs/commons/home/tlin/output/kunkle_all_2/finemap/max_snp_'
#prefix='all_anno_2'
#path='/gpfs/commons/home/tlin/output/bellenguez/bellenguez_all_2/finemap_snpvar_constrained/max_snp_'
#prefix='finemap_bellenguez_all_2'
#path='/gpfs/commons/home/tlin/output/bellenguez/bellenguez_bl/finemap_susie/max_snp_'
#path='/gpfs/commons/home/tlin/output/bellenguez/bellenguez_all_2/finemap_snpvar_constrained/max_snp_'
