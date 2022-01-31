import pandas as pd
data = pd.read_csv("/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/bellenguez_2021/bellenguez_2021_final.tsv.gz", sep = '\t', compression="gzip")

data["N"] = data.N_cases + data.N_controls
data.rename(columns={"ALLELE1": "A1", "ALLELE2": "A2","RSID":"SNP"})
data['CHR'] = data['CHR'].str.replace('chr', '')
data.to_csv("~/data/bellenguez_2021_final.tsv.gz", index = False,compression="gzip", sep = '\t' ) 
