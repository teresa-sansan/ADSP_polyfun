import pandas as pd
#print("Parquet file input:")
file='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/analysis_05_05_2024/original_sumstats/flipped/Bellenguez_et_al_2021_stage1_hg38_flipped_munged_MAF_v2.parquet'
#file = input()
print('reading file...')
df = pd.read_parquet(file)
df.to_csv('/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/ADSP_reference_panel/fine_mapping/annotations_dl/teresa/bellenguez_hg38_parquet_flipped.tsv', index = False, sep = '\t')
print(df.shape)
print('done')

