import pandas as pd
import numpy as np

polypred_max1=pd.read_csv("/gpfs/commons/home/tlin/output/bellenguez/bellenguez_all_2/polypred/adjust_beta/polypred_adj_beta_bellenguez_all2_1.prs",sep='\t')
polypred_max3=pd.read_csv("/gpfs/commons/home/tlin/output/bellenguez/bellenguez_all_2/polypred/adjust_beta/polypred_adj_beta_bellenguez_all2_3.prs",sep='\t')
polypred_max5=pd.read_csv("/gpfs/commons/home/tlin/output/bellenguez/bellenguez_all_2/polypred/adjust_beta/polypred_adj_beta_bellenguez_all2_5.prs",sep='\t')
polypred_max7=pd.read_csv("/gpfs/commons/home/tlin/output/bellenguez/bellenguez_all_2/polypred/adjust_beta/polypred_adj_beta_bellenguez_all2_7.prs",sep='\t')
polypred_max10=pd.read_csv("/gpfs/commons/home/tlin/output/bellenguez/bellenguez_all_2/polypred/adjust_beta/polypred_adj_beta_bellenguez_all2_10.prs",sep='\t')


phenotype=pd.read_csv("/gpfs/commons/home/tlin/data/phenotype_0219.2021.tsv", sep = '\t')

def combine_phenotype_polypred(polypred, phenotype = phenotype):
    input_shape = polypred.shape
    polypred=polypred.merge(phenotype, left_on = "IID", right_on = "SampleID").drop(["IID"],axis=1)
    polypred=polypred.loc[:,["SUBJID","SampleID","FID","PRS","Diagnosis","Sex","Age","APOE","Race"]]
    print("before merging:%s, after merging:%s"%(input_shape[0],polypred.shape[0]))

    return(polypred)




print("check match/unmatch:")

combine_phenotype_polypred(polypred_max1).to_csv("/gpfs/commons/home/tlin/output/prs/bellenguez_all_2/adjbeta_prs_diagnosis_0219.2021_max_snp_1.tsv",index=False, sep = '\t')
combine_phenotype_polypred(polypred_max3).to_csv("/gpfs/commons/home/tlin/output/prs/bellenguez_all_2/adjbeta_prs_diagnosis_0219.2021_max_snp_3.tsv",index=False, sep = '\t')
combine_phenotype_polypred(polypred_max5).to_csv("/gpfs/commons/home/tlin/output/prs/bellenguez_all_2/adjbeta_prs_diagnosis_0219.2021_max_snp_5.tsv",index=False, sep = '\t')
combine_phenotype_polypred(polypred_max7).to_csv("/gpfs/commons/home/tlin/output/prs/bellenguez_all_2/adjbeta_prs_diagnosis_0219.2021_max_snp_7.tsv",index=False, sep = '\t')
combine_phenotype_polypred(polypred_max10).to_csv("/gpfs/commons/home/tlin/output/prs/bellenguez_all_2/adjbeta_prs_diagnosis_0219.2021_max_snp_10.tsv",index=False, sep = '\t')

