cd /gpfs/commons/home/tlin/polyfun_omer_repo

python aggregate_finemapper_results_modified.py \
	--out-prefix /gpfs/commons/home/tlin/output/kunkle_all/finemap/all_anno \
	--sumstats /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/kunkle_2019/AD_Kunkle_etal_Stage1.parquet \
	--out ~/output/kunkle_all/finemap/aggregrate.all.txt.gz \
	--allow-missing
