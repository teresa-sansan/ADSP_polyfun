import pandas as pd

prs_001 = pd.read_csv("/gpfs/commons/home/tlin/data/plink_tutorial/height_test/change_base_file/height_pT_nospace.0.001.profile", sep = ' ').loc[:,["IID","SCORE"]]
prs_005 = pd.read_csv("/gpfs/commons/home/tlin/data/plink_tutorial/height_test/change_base_file/height_pT_nospace.0.005.profile", sep = ' ').loc[:,["IID","SCORE"]]
prs_01 = pd.read_csv("/gpfs/commons/home/tlin/data/plink_tutorial/height_test/change_base_file/height_pT_nospace.0.01.profile", sep = ' ').loc[:,["IID","SCORE"]]
prs_05 = pd.read_csv("/gpfs/commons/home/tlin/data/plink_tutorial/height_test/change_base_file/height_pT_nospace.0.05.profile", sep = ' ').loc[:,["IID","SCORE"]]
prs_1 = pd.read_csv("/gpfs/commons/home/tlin/data/plink_tutorial/height_test/change_base_file/height_pT_nospace.0.1.profile", sep = ' ').loc[:,["IID","SCORE"]]
prs_5 = pd.read_csv("/gpfs/commons/home/tlin/data/plink_tutorial/height_test/change_base_file/height_pT_nospace.0.5.profile", sep = ' ').loc[:,["IID","SCORE"]]
prs_e8 = pd.read_csv("/gpfs/commons/home/tlin/data/plink_tutorial/height_test/change_base_file/height_pT_nospace.e-8.profile", sep = ' ').loc[:,["IID","SCORE"]]
prs_e5 = pd.read_csv("/gpfs/commons/home/tlin/data/plink_tutorial/height_test/change_base_file/height_pT_nospace.e-5.profile", sep = ' ').loc[:,["IID","SCORE"]]

pheno = pd.read_csv("/gpfs/commons/home/tlin/output/prs/bellenguez_all_2/with_PC/UPDATEprs_diagnosis_0219.2021_max_snp_10_subset.tsv", sep='\t')
prs = [prs_e8, prs_e5, prs_001, prs_005, prs_01, prs_05, prs_1, prs_5]  
from functools import reduce  

prs_merge = reduce(lambda left, right:pd.merge(left,right,on=["IID"]),prs)
prs_merge = prs_merge.set_axis( ["IID","PRS_e8","PRS_e5","PRS_001","PRS_005","PRS_01","PRS_05","PRS_1","PRS_5"], axis='columns')
prs_merge.to_csv("/gpfs/commons/home/tlin/data/plink_tutorial/height_test/change_base_file/height_prs.tsv",index = False, sep='\t')


all_merge = pd.merge(pheno,prs_merge, left_on ="SampleID", right_on='IID')  
all_merge.to_csv("/gpfs/commons/home/tlin/data/plink_tutorial/height_test/change_base_file/height_pheno.tsv",index = False, sep='\t')
## cant be merged because there is no overlap
