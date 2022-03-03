import pandas as pd

path='/gpfs/commons/home/tlin/output/kunkle/kunkle_qc'
save_name='kunkle/fixed_0224/APOE'


no2SNP = pd.read_csv("/gpfs/commons/home/tlin/output/kunkle/kunkle_qc/updated_0224_remove_2SNP.prs", sep = ' ', names = ["IID","PRS"]) 
no2SNP_qc = pd.read_csv("/gpfs/commons/home/tlin/output/kunkle/kunkle_qc/updated_0224_remove_2SNP_noqc.prs", sep = ' ', names = ["IID","PRS"]) 
only2SNP = pd.read_csv("/gpfs/commons/home/tlin/output/kunkle/kunkle_qc/updated_0224_2SNP_qc.prs", sep = ' ', names = ["IID","PRS"])


pheno = pd.read_csv("/gpfs/commons/home/tlin/output/prs/bellenguez_all_2/with_PC/UPDATEprs_diagnosis_0219.2021_max_snp_10_subset.tsv", sep='\t')
prs = [only2SNP, no2SNP, no2SNP_qc]  
from functools import reduce  
prs_merge = reduce(lambda left, right:pd.merge(left,right,on=["IID"]),prs)

prs_merge = prs_merge.set_axis( ["IID","only2SNP","no2SNP","no2SNP_qc"], axis='columns')
all_merge = pd.merge(pheno,prs_merge, left_on ="SampleID", right_on='IID')
save="/gpfs/commons/home/tlin/output/prs/"+ save_name+'.tsv'
all_merge.to_csv(save,index = False, sep='\t')
