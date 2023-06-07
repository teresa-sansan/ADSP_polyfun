import pandas as pd
from functools import reduce  

#path='/gpfs/commons/home/tlin/output/wightman/prscs/'
path='/gpfs/commons/home/tlin/output/wightman/prscs/all_anno/PIP_not0/'
#path='/gpfs/commons/home/tlin/output/wightman/prscs/all_except_enformer/'
#path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_preview/PRS_hg38/prscs/plink_output/'
name='prscs'
prs_e5 = pd.read_csv(path+name+"_e-5_prs.tsv", sep = ' ', usecols = ["IID","SCORE"])
prs_001 = pd.read_csv(path+name+"_0.001_prs.tsv", sep = ' ', usecols = ["IID","SCORE"])
prs_005 = pd.read_csv(path+name+"_0.005_prs.tsv", sep = ' ', usecols = ["IID","SCORE"])   
prs_01 = pd.read_csv(path+name+"_0.01_prs.tsv", sep = ' ', usecols = ["IID","SCORE"]) 
prs_05 = pd.read_csv(path+name+"_0.05_prs.tsv", sep = ' ', usecols= ["IID","SCORE"])
prs_1 = pd.read_csv(path+name+"_0.1_prs.tsv", sep = ' ', usecols = ["IID","SCORE"]) 
prs_5 = pd.read_csv(path+name+"_0.5_prs.tsv", sep = ' ', usecols = ["IID","SCORE"]) 

# prs_e5 = pd.read_csv(path+name+"_e-5.tsv", sep = ' ', names = ["IID","PRS"])
# prs_001 = pd.read_csv(path+name+"_0.001.tsv", sep = ' ', names = ["IID","PRS"])
# prs_005 = pd.read_csv(path+name+"_0.005.tsv", sep = ' ', names = ["IID","PRS"])   
# prs_01 = pd.read_csv(path+name+"_0.01.tsv", sep = ' ', names = ["IID","PRS"]) 
# prs_05 = pd.read_csv(path+name+"_0.05.tsv", sep = ' ', names = ["IID","PRS"])
# prs_1 = pd.read_csv(path+name+"_0.1.tsv", sep = ' ', names = ["IID","PRS"]) 
# prs_5 = pd.read_csv(path+name+"_0.5.tsv", sep = ' ', names = ["IID","PRS"]) 



pheno = pd.read_csv("/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/phenotype_data_10_28_2021/all_phenotypes_unique_ancestry_subset.tsv", sep='\t')
#pheno = pd.read_csv('/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/phenotype_file/release_36K/pheno_fin.tsv', sep = '\t')

prs = [prs_e5, prs_001, prs_005, prs_01, prs_05, prs_1, prs_5]  
prs_merge = reduce(lambda left, right:pd.merge(left,right,on=["IID"]),prs)
prs_merge = prs_merge.set_axis( ["IID","PRS_e5","PRS_001","PRS_005","PRS_01","PRS_05","PRS_1","PRS_5"], axis='columns')


all_merge = pd.merge(pheno,prs_merge, left_on ="SampleID", right_on='IID')
all_merge.fillna('-1', inplace=True) # Race and Ethniciity

# save=path+name+'_36k.tsv'
save = path+name+'_pipNOT0.tsv'
all_merge.to_csv(save,index = False, sep='\t')
print("save prs to " + save)
