bl='/gpfs/commons/groups/knowles_lab/data/ldsc/polyfun/baselineLF2.2.UKB'
sumstat='/gfps/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/'

python polyfun.py \
  --compute-h2-L2 \
  --output-prefix output/jansenetal/bl_jansen \
  --sumstats $sumstat/jansenetal_2019/AD_sumstats_Jansenetal_2019sept.parquet \ 
  --ref-ld-chr $bl/baselineLF2.2.UKB. \
  --w-ld-chr $bl/weights.UKB. \
  --allow-missing


#python polyfun.py --compute-h2-L2 --output-prefix output/jansenetal/bl_jansen --sumstats /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/jansenetal_2019/AD_sumstats_Jansenetal_2019sept.parquet --ref-ld-chr /gpfs/commons/groups/knowles_lab/data/ldsc/polyfun/baselineLF2.2.UKB/baselineLF2.2.UKB. --w-ld-chr /gpfs/commons/groups/knowles_lab/data/ldsc/polyfun/baselineLF2.2.UKB/weights.UKB. --allow-missing
