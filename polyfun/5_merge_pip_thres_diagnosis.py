import pandas as pd
from functools import reduce  
import os

path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/PRS/36k_hg38/polyfun/pip_thres/middle_file'

## susie doesht have p < 0.0006
for anno in ["susie","baseline","omics","omics_dl"]:
    prs_0 = pd.read_csv(path+ '/' +anno+"_pip_0.prs", sep = ' ', names = ["IID","PRS"])
    prs_1 = pd.read_csv(path+ '/' +anno+"_pip_0.1.prs", sep = ' ', names = ["IID","PRS"])
    prs_2 = pd.read_csv(path+ '/' +anno+"_pip_0.2.prs", sep = ' ', names = ["IID","PRS"])
    prs_3 = pd.read_csv(path+ '/' +anno+"_pip_0.3.prs", sep = ' ', names = ["IID","PRS"])
    prs_4 = pd.read_csv(path+ '/' +anno+"_pip_0.4.prs", sep = ' ', names = ["IID","PRS"])
    prs_5 = pd.read_csv(path+ '/' +anno+"_pip_0.5.prs", sep = ' ', names = ["IID","PRS"])
    prs_6 = pd.read_csv(path+ '/' +anno+"_pip_0.6.prs", sep = ' ', names = ["IID","PRS"])
    prs_7 = pd.read_csv(path+ '/' +anno+"_pip_0.7.prs", sep = ' ', names = ["IID","PRS"])
    prs_8 = pd.read_csv(path+ '/' +anno+"_pip_0.8.prs", sep = ' ', names = ["IID","PRS"])
    prs_9 = pd.read_csv(path+ '/' +anno+"_pip_0.9.prs", sep = ' ', names = ["IID","PRS"])
    
    pheno=pd.read_csv('/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/phenotype_file/release_36K/pheno_ADSP_IBD.tsv', sep='\t')
    prs = [prs_0,prs_1,prs_2,prs_3,prs_4,prs_5,prs_6,prs_7,prs_8,prs_9]  
    prs_merge = reduce(lambda left, right:pd.merge(left,right,on=["IID"]),prs)
    prs_merge = prs_merge.set_axis( ["IID","PRS_0","PRS_1","PRS_2","PRS_3","PRS_4","PRS_5","PRS_6","PRS_7","PRS_8","PRS_9"], axis='columns')

    
    # prs_e4 = pd.read_csv(path+anno+"_e-4prs", sep = ' ', names = ["IID","PRS"])
    # prs_e5 = pd.read_csv(path+anno+"_e-5prs", sep = ' ', names = ["IID","PRS"])
    # #prs_e6 = pd.read_csv(path+anno+"_e-6prs", sep = ' ', names = ["IID","PRS"])
    # prs_001 = pd.read_csv(path+anno+"_0.001prs", sep = ' ', names = ["IID","PRS"])
    # prs_01 = pd.read_csv(path+anno+"_0.01prs", sep = ' ', names = ["IID","PRS"]) 
    # prs_1 = pd.read_csv(path+anno+"_0.1prs", sep = ' ', names = ["IID","PRS"]) 
    # pheno=pd.read_csv('/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/phenotype_file/release_36K/pheno_ADSP_IBD.tsv', sep='\t')
    # prs = [prs_e5, prs_e4, prs_001, prs_01,  prs_1]  
    # prs_merge = reduce(lambda left, right:pd.merge(left,right,on=["IID"]),prs)
    # prs_merge = prs_merge.set_axis( ["IID","PRS_e5","PRS_e4","PRS_001","PRS_01","PRS_1"], axis='columns')
    #prs_merge = prs_merge.set_axis( ["IID","PRS_e5","PRS_001","PRS_005","PRS_01","PRS_05","PRS_1","PRS_5"], axis='columns')
    all_merge = pd.merge(pheno,prs_merge, left_on ="SampleID", right_on='IID')
    all_merge.fillna('-1', inplace=True) # Race and Ethniciity
    all_merge.to_csv(os.path.dirname(path)  +anno+'.tsv',index = False, sep='\t')
    #print(path+anno+'.tsv')
