import pandas as pd

bellenguez = pd.read_parquet("/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/bellenguez_2021/Bellenguez_2021_stage1.parquet")
bellenguez_19 = pd.read_parquet("/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/kunkle_2019/bellenguez_2021_stage1_hg19.parquet") 

bellenguez["POS"] = bellenguez['CHR'].astype(str) +'_'+ bellenguez["BP"].astype(str)+bellenguez["A1"]+ '_' +bellenguez["A2"]
bellenguez_19["POS"] = bellenguez_19['CHR'].astype(str) +'_'+ bellenguez_19["BP"].astype(str) +bellenguez_19["Allele1"]+'_'+bellenguez_19["Allele2"]
bellenguez["POS"] = bellenguez["POS"].astype("string") 
bellenguez_19["POS"] = bellenguez_19["POS"].astype("string")


## required columns for .ma file: SNP, A1, A2, freq, b, se, p, N

new_bellenguez = bellenguez.merge(bellenguez_19.loc[:,["PVALUE","Effect","StdErr",'N_cases',"POS"]], on = "POS", how = 'inner' )
