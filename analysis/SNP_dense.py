import pandas as pd

kunkle = pd.read_parquet("/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/kunkle_2019/AD_Kunkle_etal_Stage1.parquet")
bellenguez = pd.read_parquet("/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/bellenguez_2021/Bellenguez_2021_stage1.parquet")

kunkle['loc'] = kunkle[['CHR','BP']].astype(str).apply(lambda x:'_'.join(x),axis=1)  
bellenguez['loc'] = bellenguez[['CHR','BP']].astype(str).apply(lambda x:'_'.join(x),axis=1)  


kunkle_common=kunkle[kunkle['loc'].isin(bellenguez['loc'])]
bellenguez_common=bellenguez[bellenguez["loc"].isin(kunkle['loc'])]

kunkle_common.shape ##(9734657, 8)
bellenguez_common.shape ##(9727307, 9)


kunkle_common.to_csv('~/data/kunkle_common.tsv',sep="\t",index=False)
bellenguez_common.to_csv('~/data/bellenguez_common.tsv',sep="\t",index=False)


common = pd.merge(kunkle,bellenguez, on = "loc")
new = common.drop(common.iloc[:,7:11].columns,axis=1)
new.columns=["CHR","BP","SNP","A1(kunkle)","A2(kunkle)",'N(kunkle)','Z(kunkle)',"A1(bellenguez)","A2(bellenguez)","MAF","N(bellenguez)",'Z(bellenguez)']
new.shape ##(9734657, 12)
common.to_csv('~/data/common_kunkle_bellenguez.tsv',sep='\t',index=False)


rigorous=pd.read_csv("/gpfs/commons/home/tlin/output/bl/finemap/max_causal_1/kunkle_rigorous.tsv", sep='\t')

rigorous[rigorous.iloc[:,1].isin(common.iloc[:,2])]
common_rigorous=common[common.iloc[:,2].isin(rigorous.iloc[:,1])]  
common_rigorous.to_csv("~/data/common_rigorous.tsv",sep='\t',index=False)


