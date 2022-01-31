import pandas as pd
bellenguez=pd.read_csv("/gpfs/commons/home/tlin/data/bellenguez_2021_final.tsv.gz", compression='gzip', sep='\t')
bellenguez_19 = pd.read_parquet("/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/bellenguez_2021/bellenguez_2021_stage1_hg19.parquet") 

#bellenguez["POS"] ='chr'+ bellenguez['CHR'].astype(str) +'_'+ bellenguez["BP"].astype(str) + '_' + bellenguez["A1"]+ '_' +bellenguez["A2"]
#bellenguez_19["POS"] = bellenguez_19['CHR'].astype(str) +'_'+ bellenguez_19["BP_37"].astype(str) + '_' + bellenguez_19["Allele1"]+'_'+bellenguez_19["Allele2"]
#bellenguez["POS"] = bellenguez["POS"].astype("string") 
#bellenguez_19["POS"] = bellenguez_19["POS"].astype("string")


## required columns for .ma file: SNP, A1, A2, freq, b, se, p, N

#new_bellenguez = bellenguez.merge(bellenguez_19.loc[:,["PVALUE","Effect","StdErr",'N_cases',"POS"]], on = "POS", how = 'inner' )
#duplicated_POS = new_bellenguez[new_bellenguez["POS"].duplicated()]["POS"]
#duplicated_rows = new_bellenguez[new_bellenguez["POS"].isin(duplicated_POS)]   

#duplicated_rows.to_csv("/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/bellenguez_2021/bellenguez_duplicated.tsv",sep='\t', index = False)



bellenguez_ma = bellenguez.loc[:,["RSID",'ALLELE1','ALLELE2','A1FREQ','EFFECT','SE','PVALUE','N']]
bellenguez_ma.columns = ['SNP','A1','A2','freq','b','se','p','N']  
bellenguez_ma.to_csv("~/data/bellenguez.ma", sep='\t', index=False)
