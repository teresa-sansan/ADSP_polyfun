## Author: Teresa Lin
## Dec 14, 2021

import pandas as pd
import numpy as np

polypred_max1 = pd.read_csv("/gpfs/commons/home/tlin/output/bellenguez/bellenguez_all_2/polypred/PLINKupdate_max1.prs",sep='\t')
polypred_max3 = pd.read_csv("/gpfs/commons/home/tlin/output/bellenguez/bellenguez_all_2/polypred/PLINKupdate_max3.prs",sep='\t')
polypred_max5 = pd.read_csv("/gpfs/commons/home/tlin/output/bellenguez/bellenguez_all_2/polypred/PLINKupdate_max5.prs",sep='\t')
polypred_max7 = pd.read_csv("/gpfs/commons/home/tlin/output/bellenguez/bellenguez_all_2/polypred/PLINKupdate_max7.prs",sep='\t')
polypred_max10 = pd.read_csv("/gpfs/commons/home/tlin/output/bellenguez/bellenguez_all_2/polypred/PLINKupdate_max10.prs",sep='\t')
#phenotype_chirag = pd.read_csv("/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/phenotype_data_10_28_2021/all_phenotypes_unique_ancestry_all_all.tsv", sep='\t')
phenotype_chirag_noNA = pd.read_csv("/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/phenotype_data_10_28_2021/all_phenotypes_unique_ancestry_subset.tsv", sep = '\t')


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


max1 = combine_phenotype_polypred(polypred_max1,"/gpfs/commons/home/tlin/output/prs/bellenguez_all_2/with_PC/UPDATEprs_diagnosis_0219.2021_max_snp_1_subset.tsv")
max3 = combine_phenotype_polypred(polypred_max3,"/gpfs/commons/home/tlin/output/prs/bellenguez_all_2/with_PC/UPDATEprs_diagnosis_0219.2021_max_snp_3_subset.tsv")
max5 = combine_phenotype_polypred(polypred_max5,"/gpfs/commons/home/tlin/output/prs/bellenguez_all_2/with_PC/UPDATEprs_diagnosis_0219.2021_max_snp_5_subset.tsv")
max7 = combine_phenotype_polypred(polypred_max7,"/gpfs/commons/home/tlin/output/prs/bellenguez_all_2/with_PC/UPDATEprs_diagnosis_0219.2021_max_snp_7_subset.tsv")
max10 = combine_phenotype_polypred(polypred_max10,"/gpfs/commons/home/tlin/output/prs/bellenguez_all_2/with_PC/UPDATEprs_diagnosis_0219.2021_max_snp_10_subset.tsv")
