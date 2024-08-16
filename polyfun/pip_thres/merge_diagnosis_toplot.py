import pandas as pd
from functools import reduce  
import sys

anno=sys.argv[1]
path=sys.argv[2]
flag=sys.argv[3]

if (flag='no_qc'):
    print('geno no qc')
    path=path+flag


pheno=pd.read_csv('/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/phenotype_file/release_36K/pheno_ADSP_IBD.tsv', sep='\t')

prs_0=pd.read_csv(path+"middle_file/%s_pip_0.prs"%anno ,sep = ' ', names = ["IID","PRS"])
prs_1=pd.read_csv(path+"middle_file/%s_pip_0.1.prs"%anno ,sep = ' ', names = ["IID","PRS"])
prs_2=pd.read_csv(path+"middle_file/%s_pip_0.2.prs"%anno ,sep = ' ', names = ["IID","PRS"])
prs_3=pd.read_csv(path+"middle_file/%s_pip_0.3.prs"%anno ,sep = ' ', names = ["IID","PRS"])
prs_4=pd.read_csv(path+"middle_file/%s_pip_0.4.prs"%anno ,sep = ' ', names = ["IID","PRS"])
prs_5=pd.read_csv(path+"middle_file/%s_pip_0.5.prs"%anno ,sep = ' ', names = ["IID","PRS"])
prs_6=pd.read_csv(path+"middle_file/%s_pip_0.6.prs"%anno ,sep = ' ', names = ["IID","PRS"])
prs_7=pd.read_csv(path+"middle_file/%s_pip_0.7.prs"%anno ,sep = ' ', names = ["IID","PRS"])
prs_8=pd.read_csv(path+"middle_file/%s_pip_0.8.prs"%anno ,sep = ' ', names = ["IID","PRS"])
prs_9=pd.read_csv(path+"middle_file/%s_pip_0.9.prs"%anno ,sep = ' ', names = ["IID","PRS"])

prs = [prs_0, prs_1, prs_2, prs_3, prs_4, prs_5, prs_6, prs_7, prs_8, prs_9]        
prs_merge = reduce(lambda left, right:pd.merge(left,right,on=["IID"]),prs)
prs_merge = prs_merge.set_axis( ["IID","PRS_0","PRS_1","PRS_2","PRS_3","PRS_4","PRS_5","PRS_6","PRS_7","PRS_8","PRS_9"], axis='columns')
all_merge = pd.merge(pheno,prs_merge, left_on ="SampleID", right_on='IID')
all_merge.fillna('-1', inplace=True) # Race and Ethniciity
save= path + 'prs_%s.tsv'%anno


all_merge.to_csv(save,index = False, sep='\t')
print("save prs to " + save)
