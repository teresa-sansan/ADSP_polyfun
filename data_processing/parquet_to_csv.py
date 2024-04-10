import pandas as pd
import os,sys

#file = sys.argv[1]
#file = '/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/bellenguez_2021/AD_Bellenguez_polyfun_fixed.parquet'
# file = '/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/wightman_2021/wightman_2021_fixed.parquet'
# file = '/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/kunkle_2019/AD_Kunkle_etal_Stage1.parquet'
file = '/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/jansenetal_2019/AD_sumstats_Jansenetal_2019sept.parquet'

output_path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/sumstat_chirag/'
df = pd.read_parquet(file)
df.to_csv(output_path + 'jansen' +'.tsv', sep = '\t', index = False)
print('write %s into tsv'%file)
print(output_path + 'bellenguez')