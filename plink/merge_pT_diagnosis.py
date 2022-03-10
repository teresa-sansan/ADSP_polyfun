import pandas as pd

path='/gpfs/commons/home/tlin/output/cT/kunkle/fixed_0224/qc_on_base/'
#path='/gpfs/commons/home/tlin/output/cT/kunkle/qc_check/'
save_name='kunkle/fixed_0224/kunkle_qc_on_base'

prs_e5 = pd.read_csv(path+"pT_e-5.prs", sep = ' ', names = ["IID","PRS"])
prs_001 = pd.read_csv(path+"pT_0.001.prs", sep = ' ', names = ["IID","PRS"])
prs_005 = pd.read_csv(path+"pT_0.005.prs", sep = ' ', names = ["IID","PRS"])   
prs_01 = pd.read_csv(path+"pT_0.01.prs", sep = ' ', names = ["IID","PRS"]) 
prs_05 = pd.read_csv(path+"pT_0.05.prs", sep = ' ', names = ["IID","PRS"])
prs_1 = pd.read_csv(path+"pT_0.1.prs", sep = ' ', names = ["IID","PRS"]) 
prs_5 = pd.read_csv(path+"pT_0.5.prs", sep = ' ', names = ["IID","PRS"]) 


pheno = pd.read_csv("/gpfs/commons/home/tlin/output/prs/bellenguez_all_2/with_PC/UPDATEprs_diagnosis_0219.2021_max_snp_10_subset.tsv", sep='\t')
prs = [prs_e5, prs_001, prs_005, prs_01, prs_05, prs_1, prs_5]  
from functools import reduce  
prs_merge = reduce(lambda left, right:pd.merge(left,right,on=["IID"]),prs)

prs_merge = prs_merge.set_axis( ["IID","PRS_e5","PRS_001","PRS_005","PRS_01","PRS_05","PRS_1","PRS_5"], axis='columns')
all_merge = pd.merge(pheno,prs_merge, left_on ="SampleID", right_on='IID')
save="/gpfs/commons/home/tlin/output/prs/"+ save_name+'.tsv'
all_merge.to_csv(save,index = False, sep='\t')
