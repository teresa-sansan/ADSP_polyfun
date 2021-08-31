cd /gpfs/commons/home/tlin/polyfun_omer_repo  

prefix='/gpfs/commons/home/tlin/output/kunkle_all/finemap/all_anno.22'

for i in prefix
do
echo $i
done

python aggregate_finemapper_results.py \
	--out-prefix $prefix \
	--sumstat /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/kunkle_2019/AD_Kunkle_etal_Stage1.parquet \
	--out /gpfs/commons/home/tlin/output/kunkle_all/finemap/kunkle_agg.txt.gz \
	--chr 22 \
	--allow-missing-job \
 	--adjust-beta-freq
