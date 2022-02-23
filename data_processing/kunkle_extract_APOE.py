import pandas as pd
kunkle = pd.read_parquet('/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/kunkle_2019/AD_Kunkle_etal_Stage1.parquet')
APOE = kunkle[(kunkle.CHR == 19) &(kunkle.BP>45409011)&(kunkle.BP<45412650)] 
APOE.to_csv("/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/kunkle_2019/APOE_only_no_qc.tsv", sep = '\t', index = False)
