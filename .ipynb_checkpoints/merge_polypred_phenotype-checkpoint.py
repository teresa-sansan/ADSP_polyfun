import pandas as pd
import numpy as np

#polypred_max1 = pd.read_csv("/gpfs/commons/home/tlin/output/bellenguez/bellenguez_all_2/polypred/PLINKupdate_max1.prs",sep='\t')
#polypred_max3 = pd.read_csv("/gpfs/commons/home/tlin/output/bellenguez/bellenguez_all_2/polypred/PLINKupdate_max3.prs",sep='\t')
#polypred_max5 = pd.read_csv("/gpfs/commons/home/tlin/output/bellenguez/bellenguez_all_2/polypred/PLINKupdate_max5.prs",sep='\t')
#polypred_max7 = pd.read_csv("/gpfs/commons/home/tlin/output/bellenguez/bellenguez_all_2/polypred/PLINKupdate_max7.prs",sep='\t')
polypred_max10 = pd.read_csv("/gpfs/commons/home/tlin/output/bellenguez/bellenguez_all_2/polypred/PLINKupdate_max10.prs",sep='\t')
phenotype_chirag = pd.read_csv("/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/phenotype_data_10_28_2021/all_phenotypes_unique_ancestry_all_all.tsv", sep='\t')
phenotype_chirag_noNA = pd.read_csv("/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/phenotype_data_10_28_2021/all_phenotypes_unique_ancestry_subset.tsv", sep = '\t')


print(phenotype_chirag_noNA.shape)


def combine_phenotype_polypred(polypred, phenotype = phenotype_chirag_noNA):
    input_shape = phenotype.shape
    polypred=phenotype.merge(polypred, right_on = "IID", left_on = "SampleID").drop(["IID"],axis=1)
    polypred=polypred.loc[:,["SUBJID","SampleID","PRS","AD_status_final","Sex","age_covariate","APOE","Race","Ethnicity","final_population","X1","X2","X3","X4","X5","X6","X7","X8","X9","X10","X11","X12","X13","X14","X15","X16","X17","X18","X19","X20"]]
    polypred = polypred.rename(columns={"AD_status_final":"Diagnosis","age_covariate":"Age"}) 
    print("before merging:%s, after merging:%s"%(input_shape[0],polypred.shape[0]))

    return(polypred)


# def combine_phenotype_polypred(polypred, phenotype = phenotype_chirag_noNA):
#     input_shape = polypred.shape
#     polypred=polypred.merge(phenotype, left_on = "IID", right_on = "SampleID").drop(["IID"],axis=1)
#     polypred=polypred.loc[:,["SUBJID","SampleID","PRS","AD_status_final","Sex","age_covariate","APOE","Race","Ethnicity","final_population","X1","X2","X3","X4","X5","X6","X7","X8","X9","X10","X11","X12","X13","X14","X15","X16","X17","X18","X19","X20"]]
#     polypred = polypred.rename(columns={"AD_status_final":"Diagnosis","age_covariate":"Age"}) 
#     print("before merging:%s, after merging:%s"%(input_shape[0],polypred.shape[0]))

#     return(polypred)






#max1 = combine_phenotype_polypred(polypred_max1)
#max3 = combine_phenotype_polypred(polypred_max3)
#max5 = combine_phenotype_polypred(polypred_max5)
#max7 = combine_phenotype_polypred(polypred_max7)
max10 = combine_phenotype_polypred(polypred_max10)


##because idk why but direcly merging them together will generate strange values ( they probably just randomly filled up the NAs)

max10['Race'] = max10['Race'].replace(np.nan, -100)
max10['Ethnicity'] = max10['Ethnicity'].replace(np.nan, -100)    

#max1.to_csv("/gpfs/commons/home/tlin/output/prs/bellenguez_all_2/with_PC/UPDATEprs_diagnosis_0219.2021_max_snp_1_subset.tsv",index=False, sep = '\t')
#max3.to_csv("/gpfs/commons/home/tlin/output/prs/bellenguez_all_2/with_PC/UPDATEprs_diagnosis_0219.2021_max_snp_3_subset.tsv",index=False, sep = '\t')
#max5.to_csv("/gpfs/commons/home/tlin/output/prs/bellenguez_all_2/with_PC/UPDATEprs_diagnosis_0219.2021_max_snp_5_subset.tsv",index=False, sep = '\t')
#max7.to_csv("/gpfs/commons/home/tlin/output/prs/bellenguez_all_2/with_PC/UPDATEprs_diagnosis_0219.2021_max_snp_7_subset.tsv",index=False, sep = '\t')
max10.to_csv("/gpfs/commons/home/tlin/output/prs/bellenguez_all_2/with_PC/UPDATEprs_diagnosis_0219.2021_max_snp_10_subset_new.tsv",index=False, sep = '\t')





