## Author: Teresa Lin
## Dec 14, 2021

import pandas as pd
import numpy as np
from functools import reduce
#save_name='polypred/wightman/fixed_0224'
#path='/gpfs/commons/home/tlin/output/wightman/fixed_0224/finemap/polypred/'
path='/gpfs/commons/home/tlin/output/bellenguez/bellenguez_fixed_0224/finemap/try_rescue_not_converge/polypred/polypred_genomewide.tsv.prs'
#polypred_max1 = pd.read_csv(path+'max_snp_1.wightman_polypred.tsv.prs',sep='\t',names=["FID","IID","PRS"] )
#polypred_max3 = pd.read_csv(path+'max_snp_3.wightman_polypred.tsv.prs',sep='\t',names=["FID","IID","PRS"])
#polypred_max5 = pd.read_csv(path+'max_snp_5.wightman_polypred.tsv.prs',sep='\t',names=["FID","IID","PRS"])
#polypred_max7 = pd.read_csv(path+'max_snp_7.wightman_polypred.tsv.prs',sep='\t',names=["FID","IID","PRS"])
polypred_max10 = pd.read_csv(path,sep='\t',names=["FID","IID","PRS"])
phenotype_chirag = pd.read_csv("/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/phenotype_data_10_28_2021/all_phenotypes_unique_ancestry_all_all.tsv", sep='\t')
phenotype_chirag_noNA = pd.read_csv("/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/phenotype_data_10_28_2021/all_phenotypes_unique_ancestry_subset.tsv", sep = '\t')

#prs=[polypred_max1,polypred_max3,polypred_max5,polypred_max7,polypred_max10]
#prs_merge = reduce(lambda left, right:pd.merge(left,right,on=["IID"]),prs)
#print(prs_merge)

#prs_merge = prs_merge.set_axis( ["IID","PRS_max_snp_1","PRS_max_snp_3","PRS_max_snp_5","PRS_max_snp_7","PRS_max_snp_10"],axis='columns')
#print(prs_merge)

#all_merge = pd.merge(pheno,prs_merge, left_on ="SampleID", right_on='IID')
#save="/gpfs/commons/home/tlin/output/prs/"+ save_name+'.tsv'
#all_merge.to_csv(save,index = False, sep='\t')


print(phenotype_chirag_noNA.shape)
#print("check if there are NAs")
#print(phenotype_chirag_noNA.isna().sum())

def combine_phenotype_polypred(polypred, savename, phenotype = phenotype_chirag_noNA):
    input_shape = phenotype.shape
    polypred=phenotype.merge(polypred, right_on = "IID", left_on = "SampleID")
    polypred = polypred.drop(columns = ['FID','IID','flag_age_covariate','Duplicate_SUBJID'],axis = 1)   ## removed those unwanted columns. 
    polypred = polypred.rename(columns={"AD_status_final":"Diagnosis","age_covariate":"Age"})  ## rename columns for simplification
    print("before merging:%s, after merging:%s"%(input_shape[0],polypred.shape[0]))
    print(polypred.isna().sum())
    ## because idk why but when leaving them with NAs will generate wrong result when importing to R.
    ## So I directly filled them up with -100. 
    polypred['Race'] = polypred['Race'].replace(np.nan, -100)
    polypred['Ethnicity'] = polypred['Ethnicity'].replace(np.nan, -100)
    polypred.to_csv(savename, index=False, sep ='\t')
    print("save to " , savename)
    return(polypred)


#max1 = combine_phenotype_polypred(polypred_max1,"/gpfs/commons/home/tlin/output/prs/bellenguez_all_2/with_PC/UPDATEprs_diagnosis_0219.2021_max_snp_1_subset.tsv")
#max3 = combine_phenotype_polypred(polypred_max3,"/gpfs/commons/home/tlin/output/prs/bellenguez_all_2/with_PC/UPDATEprs_diagnosis_0219.2021_max_snp_3_subset.tsv")
#max5 = combine_phenotype_polypred(polypred_max5,"/gpfs/commons/home/tlin/output/prs/bellenguez_all_2/with_PC/UPDATEprs_diagnosis_0219.2021_max_snp_5_subset.tsv")
#max7 = combine_phenotype_polypred(polypred_max7,"/gpfs/commons/home/tlin/output/prs/bellenguez_all_2/with_PC/UPDATEprs_diagnosis_0219.2021_max_snp_7_subset.tsv")
max10 = combine_phenotype_polypred(polypred_max10,"/gpfs/commons/home/tlin/output/prs/polypred/bellenguez/fixed_0224_polypred.prs")
