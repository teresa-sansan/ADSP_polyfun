import pandas as pd
from functools import reduce  
import sys

path=sys.argv[1]
save_name=sys.argv[2]
    
pheno=pd.read_csv('/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/phenotype_file/release_36K/pheno_ADSP_IBD.tsv', sep='\t')

prs_e4 = pd.read_csv(path + "pT_e-4.prs", sep = ' ', names = ["IID","PRS"])
prs_e5 = pd.read_csv(path + "pT_e-5.prs", sep = ' ', names = ["IID","PRS"])
prs_e6 = pd.read_csv(path + "pT_e-6.prs", sep = ' ', names = ["IID","PRS"])
prs_001 = pd.read_csv(path + "pT_0.001.prs", sep = ' ', names = ["IID","PRS"])
prs_01 = pd.read_csv(path + "pT_0.01.prs", sep = ' ', names = ["IID","PRS"]) 
prs_1 = pd.read_csv(path + "pT_0.1.prs", sep = ' ', names = ["IID","PRS"]) 
prs_0 = pd.read_csv(path + "pT_0.prs", sep = ' ', names = ["IID","PRS"]) 

pheno=pd.read_csv('/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/phenotype_file/release_36K/pheno_ADSP_IBD.tsv', sep='\t')
prs = [prs_e6, prs_e5, prs_e4, prs_001, prs_01, prs_1, prs_0]  
prs_merge = reduce(lambda left, right:pd.merge(left,right,on=["IID"]),prs)
prs_merge = prs_merge.set_axis( ["IID","PRS_e6","PRS_e5","PRS_e4","PRS_001","PRS_01","PRS_1", "PRS_0"], axis='columns')

all_merge = pd.merge(pheno,prs_merge, left_on ="SampleID", right_on='IID')
all_merge.fillna('-1', inplace=True) # Race and Ethniciity
all_merge.to_csv(save_name,index = False, sep='\t')
print('save to ' + save_name)
