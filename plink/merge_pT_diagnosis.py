import pandas as pd

#path='/gpfs/commons/home/tlin/output/cT/kunkle/fixed_0224/qc_target_maf01/'
#path='/gpfs/commons/home/tlin/output/cT/kunkle/qc_check/'
#save_name='kunkle/fixed_0224/kunkle_qc_target_maf01'

#path='/gpfs/commons/home/tlin/output/cT/wightman/before_qc/'
#save_name='kunkle/fixed_0224/qc_on_variant_maf01'
path='/gpfs/commons/home/tlin/output/cT/wightman/qc_on_individual/'
save_name='wightman/qc_on_variant_update'

prs_e5 = pd.read_csv(path+"pT_e-5_update.prs", sep = ' ', names = ["IID","PRS"])
prs_001 = pd.read_csv(path+"pT_0.001_update.prs", sep = ' ', names = ["IID","PRS"])
prs_005 = pd.read_csv(path+"pT_0.005_update.prs", sep = ' ', names = ["IID","PRS"])   
prs_01 = pd.read_csv(path+"pT_0.01_update.prs", sep = ' ', names = ["IID","PRS"]) 
prs_05 = pd.read_csv(path+"pT_0.05_update.prs", sep = ' ', names = ["IID","PRS"])
prs_1 = pd.read_csv(path+"pT_0.1_update.prs", sep = ' ', names = ["IID","PRS"]) 
prs_5 = pd.read_csv(path+"pT_0.5_update.prs", sep = ' ', names = ["IID","PRS"]) 


pheno = pd.read_csv("/gpfs/commons/home/tlin/output/prs/bellenguez_all_2/with_PC/UPDATEprs_diagnosis_0219.2021_max_snp_10_subset.tsv", sep='\t')
prs = [prs_e5, prs_001, prs_005, prs_01, prs_05, prs_1, prs_5]  
from functools import reduce  
prs_merge = reduce(lambda left, right:pd.merge(left,right,on=["IID"]),prs)

prs_merge = prs_merge.set_axis( ["IID","PRS_e5","PRS_001","PRS_005","PRS_01","PRS_05","PRS_1","PRS_5"], axis='columns')
all_merge = pd.merge(pheno,prs_merge, left_on ="SampleID", right_on='IID')
save="/gpfs/commons/home/tlin/output/prs/"+ save_name+'.tsv'
all_merge.to_csv(save,index = False, sep='\t')
