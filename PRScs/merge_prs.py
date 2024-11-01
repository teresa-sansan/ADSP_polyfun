import pandas as pd
from functools import reduce  
import sys

path = sys.argv[1]
save = sys.argv[2]

name='prs'
prs_e5 = pd.read_csv(path+name+"_e-5.tsv", sep = ' ', names = ["IID","SCORE"])
prs_e4 = pd.read_csv(path+name+"_e-4.tsv", sep = ' ', names = ["IID","SCORE"])
prs_001 = pd.read_csv(path+name+"_0.001.tsv", sep = ' ', names = ["IID","SCORE"])
#prs_005 = pd.read_csv(path+name+"_0.005.tsv", sep = ' ', usecols = ["IID","SCORE"])   
prs_01 = pd.read_csv(path+name+"_0.01.tsv", sep = ' ', names = ["IID","SCORE"]) 
#prs_05 = pd.read_csv(path+name+"_0.05.tsv", sep = ' ', names= ["IID","SCORE"])
prs_1 = pd.read_csv(path+name+"_0.1.tsv", sep = ' ', names = ["IID","SCORE"]) 
prs_9 = pd.read_csv(path+name+"_0.9.tsv", sep = ' ', names = ["IID","SCORE"]) 
#prs_5 = pd.read_csv(path+name+"_0.5.tsv", sep = ' ', names = ["IID","SCORE"]) 


## 36k w. fixed IBD
pheno = pd.read_csv('/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/phenotype_file/release_36K/pheno_ADSP_IBD.tsv', sep = '\t')
prs = [ prs_e5 ,prs_e4,prs_001, prs_01 , prs_1,prs_9]  
prs_merge = reduce(lambda left, right:pd.merge(left,right,on=["IID"]),prs)
prs_merge = prs_merge.set_axis( ["IID","PRS_e5","PRS_e4","PRS_001","PRS_01","PRS_1","PRS_9"], axis='columns')


all_merge = pd.merge(pheno,prs_merge, left_on ="SampleID", right_on='IID')
all_merge.fillna('-1', inplace=True) # Race and Ethniciity

save=path+'_'+save+ '_36k.tsv'
# save=path+'prscs_wightman_ADSP_ibd_36k.tsv'
# save = path+name+'_pipNOT0.tsv'
all_merge.to_csv(save,index = False, sep='\t')
print("save prs to " + save)
