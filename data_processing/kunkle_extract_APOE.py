import pandas as pd
kunkle = pd.read_parquet('AD_Kunkle_etal_Stage1.parquet')
APOE = kunkle[(kunkle.CHR == 19) &(kunkle.BP>45409011)&(kunkle.BP<45412079)] 
APOE.to_csv("APOE_only.tsv", sep = '\t', index = False)
