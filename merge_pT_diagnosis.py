import pandas as pd
!cd /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap/pT

prs_001 = pd.read_csv("pT_0.001.prs", sep = ' ', names = ["IID","PRS"])  
prs_01 = pd.read_csv("pT_0.01.prs", sep = ' ', names = ["IID","PRS"]) 
prs_05 = pd.read_csv("pT_0.05.prs", sep = ' ', names = ["IID","PRS"])
prs_1 = pd.read_csv("pT_0.1.prs", sep = ' ', names = ["IID","PRS"]) 
prs_2 = pd.read_csv("pT_0.2.prs", sep = ' ', names = ["IID","PRS"]) 
prs_5 = pd.read_csv("pT_0.5.prs", sep = ' ', names = ["IID","PRS"]) 

pheno = pd.read_csv("/gpfs/commons/home/tlin/output/prs/bellenguez_all_2/with_PC/UPDATEprs_diagnosis_0219.2021_max_snp_10_subset.tsv", sep='\t')
prs = [prs_001, prs_01, prs_05, prs_1, prs_2, prs_5]  
from functools import reduce  
prs_merge = reduce(lambda left, right:pd.merge(left,right,on=["IID"]),prs)

prs_merge.set_axis( ["IID","PRS_001","PRS_01","PRS_05","PRS_01","PRS_02","PRS_05"], axis='columns')

all_merge = pd.merge(pheno,prs_merge, left_on ="SampleID", right_on='IID')  
all_merge.to_csv("/gpfs/commons/home/tlin/output/prs/bellenguez_all_2/pT/pT_PRS_withPC.tsv",index = False, sep='\t')
