import pandas as pd

sumstat1 = pd.read_parquet("/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/bellenguez_2021/Bellenguez_2021_stage1.parquet")
sumstat2 = pd.read_parquet("/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/bellenguez_2021/bellenguez_2021_stage1_hg19.parquet")

print(sumstat1.shape, sumstat2.shape)

sumstat1['loc'] = sumstat1[['CHR','BP']].astype(str).agg('-'.join, axis=1).str.cat(sumstat1[["A1",'A2']],sep='-')  
sumstat2['loc'] = sumstat2[['CHR','BP']].astype(str).agg('-'.join, axis=1).str.cat(sumstat1[["A1",'A2']],sep='-')  


sumstat = sumstat1.merge(sumstat2, on ='loc', how ='inner')
print("sumstat_shape=")
print(sumstat.shape)

sumstat.rename(columns = {"CHR_x": "CHR", "BP_x": "BP"},errors="raise")
sumstat.drop(["loc","CHR_y","BP_y","Allele1","Allele2"], axis=1)
sumstat.to_csv("/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/bellenguez_2021/merge_bellenguez.tsv",sep = '\t', index = False)


ma = data[["SNP","A1","A2","MAF","Effect","StdErr","PVALUE","N"]]
ma = ma.rename(columns = {"MAF": "freq", "Effect": "b","StdErr":"se", "PVALUE":"p"})
ma.to_csv("/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/bellenguez_2021/Bellenguez_2021_stage1.ma",sep = '\t', index = False)

