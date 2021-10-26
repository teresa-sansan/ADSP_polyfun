import pandas as pd
import numpy as np

polypred_max1 = pd.read_csv("/gpfs/commons/home/tlin/output/bellenguez/bellenguez_all_2/polypred/polypred_bellenguez_max1.prs",sep='\t')
polypred_max3 = pd.read_csv("/gpfs/commons/home/tlin/output/bellenguez/bellenguez_all_2/polypred/polypred_bellenguez_max3.prs",sep='\t')
polypred_max5 = pd.read_csv("/gpfs/commons/home/tlin/output/bellenguez/bellenguez_all_2/polypred/polypred_bellenguez_max5.prs",sep='\t')
polypred_max7 = pd.read_csv("/gpfs/commons/home/tlin/output/bellenguez/bellenguez_all_2/polypred/polypred_bellenguez_max7.prs",sep='\t')
polypred_max10 = pd.read_csv("/gpfs/commons/home/tlin/output/bellenguez/bellenguez_all_2/polypred/polypred_bellenguez_max10.prs",sep='\t')
manifest = pd.read_csv("/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/phenotype_data/SampleManifest_DS_2021.02.19_ALL.txt", sep='\t') ## the one with mapping info
phenotype1 = pd.read_csv("/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/phenotype_data/ADNIPhenotypes_DS_2021.02.19_ALL.txt", sep='\t')
phenotype2 = pd.read_csv("/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/phenotype_data/ADSPCaseControlPhenotypes_DS_2021.02.19_ALL.txt", sep='\t')
phenotype3 = pd.read_csv("/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/phenotype_data/ADSPFamilyBasedPhenotypes_DS_2021.02.19_ALL.txt", sep='\t')
phenotype4 = pd.read_csv("/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/phenotype_data/PSPCBDPhenotypes_DS_2021.02.19_ALL.txt", sep='\t')


phenotype=pd.concat([phenotype1, phenotype2])

def merge_polypred_phenotype(polypred, pheno=phenotype, manifest=manifest):
    subjid = manifest[manifest.SampleID.isin(polypred.IID)].SUBJID  ## map Sample ID to Subject ID
    pheno=pheno[pheno.SUBJID.isin(subjid)]
    print("pheno shape = ",pheno.shape[0])
    ID = pheno.SUBJID.drop_duplicates() ##SubjID 
    found = manifest[manifest.SUBJID.isin(ID)].SampleID ## map Subject ID back to Sample ID
    print("are there any unmapped polypred ID?")
    print("polypred shape = ", polypred.shape, "; not-matching shape = ", polypred[~polypred.IID.isin(found)].shape) 
    print("pheno before removing NA in diagnosis= ", pheno.shape)
    pheno = pheno.dropna(subset = ["PrevAD",'IncAD'])
    print("pheno after removing NA in diagnosis = ", pheno.shape)
    pheno["Diagnosis"]= pheno[['PrevAD', 'IncAD']].max(axis=1).astype(int)

    
    reverse_map = manifest[manifest.SUBJID.isin(ID)]
    polypred.rename({"IID":"SampleID"}, axis=1, inplace=True)
    info =pd.merge(pd.merge(polypred, reverse_map[reverse_map.SampleID.isin(polypred.SampleID)][["SUBJID","SampleID"]]),pheno).loc[:,["SampleID","PRS","SUBJID","Sex","Diagnosis"]]
    return(info)


max1 = merge_polypred_phenotype(polypred_max1)
max3 = merge_polypred_phenotype(polypred_max3)
max5 = merge_polypred_phenotype(polypred_max5)
max7 = merge_polypred_phenotype(polypred_max7)
max10 = merge_polypred_phenotype(polypred_max10)


max1.to_csv("/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/phenotype_data/diagnosis/prs_diagnosis_0219.2021_max_snp_1.tsv",index=False, sep = '\t')
max3.to_csv("/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/phenotype_data/diagnosis/prs_diagnosis_0219.2021_max_snp_3.tsv",index=False, sep = '\t')
max5.to_csv("/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/phenotype_data/diagnosis/prs_diagnosis_0219.2021_max_snp_5.tsv",index=False, sep = '\t')
max7.to_csv("/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/phenotype_data/diagnosis/prs_diagnosis_0219.2021_max_snp_7.tsv",index=False, sep = '\t')
max10.to_csv("/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/phenotype_data/diagnosis/prs_diagnosis_0219.2021_max_snp10.tsv",index=False, sep = '\t')





