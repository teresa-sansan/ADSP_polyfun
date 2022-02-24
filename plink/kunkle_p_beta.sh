import pandas as pd
kunkle = pd.read_parquet("kunkle_2019_teresa/AD_Kunkle_etal_Stage1.parquet")
beta = pd.read_csv("/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/kunkle_2019/Kunkle_etal_Stage1_results_modified.txt", sep = ' ')

beta["pos"] = beta.Chromosome.astype("str")+'_'+beta.Position.astype("str")+'_'+beta.A1.astype("str")+'_'+beta.A2.astype("str") 
kunkle["pos"]=kunkle.CHR.astype("str")+'_'+kunkle.BP.astype("str")+'_'+kunle.A1.astype("str")+'_'+kunkle.A2.astype("str")

data = pd.merge(kunkle, beta.loc[:,["Beta","SE","Pvalue",'pos']], on='pos')
data = data.drop(columns = 'pos')   
data.to_csv("/gpfs/commons/home/tlin/data/kunkle_etal_Stage1_info.tsv",sep= '\t',index = False) 
