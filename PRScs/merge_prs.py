import pandas as pd
from functools import reduce  

path='/gpfs/commons/home/tlin/output/wightman/prscs/'

prs_e5 = pd.read_csv(path+"p_e-5_prs.tsv", sep = ' ', usecols = ["IID","SCORE"])
prs_001 = pd.read_csv(path+"p_0.001_prs.tsv", sep = ' ', usecols = ["IID","SCORE"])
prs_005 = pd.read_csv(path+"p_0.005_prs.tsv", sep = ' ', usecols = ["IID","SCORE"])   
prs_01 = pd.read_csv(path+"p_0.01_prs.tsv", sep = ' ', usecols = ["IID","SCORE"]) 
prs_05 = pd.read_csv(path+"p_0.05_prs.tsv", sep = ' ', usecols= ["IID","SCORE"])
prs_1 = pd.read_csv(path+"p_0.1_prs.tsv", sep = ' ', usecols = ["IID","SCORE"]) 
prs_5 = pd.read_csv(path+"p_0.5_prs.tsv", sep = ' ', usecols = ["IID","SCORE"]) 

pheno = pd.read_csv("/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/phenotype_data_10_28_2021/all_phenotypes_unique_ancestry_subset.tsv", sep='\t')
prs = [prs_e5, prs_001, prs_005, prs_01, prs_05, prs_1, prs_5]  
prs_merge = reduce(lambda left, right:pd.merge(left,right,on=["IID"]),prs)
prs_merge = prs_merge.set_axis( ["IID","PRS_e5","PRS_0.001","PRS_0.005","PRS_0.01","PRS_0.05","PRS_0.1","PRS_0.5"], axis='columns')

all_merge = pd.merge(pheno,prs_merge, left_on ="SampleID", right_on='IID')
all_merge.fillna('-1', inplace=True) # Race and Ethniciity
save=path+'prscs_17k.tsv'
all_merge.to_csv(save,index = False, sep='\t')
print("save prs to " + save)