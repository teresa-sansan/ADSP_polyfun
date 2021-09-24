cd /gpfs/commons/home/tlin/polyfun_omer_repo

kunkle='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/kunkle_2019/AD_Kunkle_etal_Stage1.parquet'
bellenguez='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/bellenguez_2021/Bellenguez_2021_stage1.parquet'
path='/gpfs/commons/home/tlin/output/kunkle_all_2/finemap/max_snp_1'
prefix='all_anno'

for chr in {1..22}
do
	python aggregate_finemapper_results_modified.py \
		--out-prefix $path/$prefix \
		--sumstats $kunkle \
		--out $path/chr${chr}.aggregrate.all.txt.gz \
		--allow-missing \
       		--chr $chr
done


#prefix='/gpfs/commons/home/tlin/output/kunkle_all/finemap/all_anno'
#output='~/output/kunkle_all/finemap/aggregrate.all.txt.gz'
#prefix='/gpfs/commons/home/tlin/output/bellenguez/bellenguez_all_2/finemap/max_snp_1/finemap_bellenguez_all_2'
