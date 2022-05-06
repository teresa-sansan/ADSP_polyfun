import pandas as pd

#path = '/gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap/polypred/'
#path = '/gpfs/commons/home/tlin/output/bellenguez/bellenguez_all_2/polypred/PLINKupdate_max'
#path = '/gpfs/commons/home/tlin/output/prs/bellenguez_susie/susie_prs_diagnosis_0219.2021_max_snp_'
#path = '/gpfs/commons/home/tlin/output/prs/bellenguez_all_2/with_PC/UPDATEprs_diagnosis_0219.2021_max_snp_'
#path = '/gpfs/commons/home/tlin/output/bellenguez/bellenguez_fixed_0224/finemap/polypred/'
#path = '/gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap/polypred/'
path = '/gpfs/commons/home/tlin/output/kunkle/kunkle_fixed_0224/finemap/polypred/'
save = '/gpfs/commons/home/tlin/output/prs/kunkle/fixed_0224/'


prs1 = pd.read_csv(path+'1.prs', sep = '\t')
prs3 = pd.read_csv(path+'3.prs', sep = '\t') 
prs5 = pd.read_csv(path+'5.prs', sep = '\t')
prs7 = pd.read_csv(path+'7.prs', sep = '\t')
prs10 = pd.read_csv(path+'10.prs', sep = '\t')


pheno = pd.read_csv('/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/phenotype_data_10_28_2021/all_phenotypes_unique_ancestry_subset.tsv', sep = '\t')
prs = pd.DataFrame({'PRS1':prs1.PRS,'PRS3':prs3.PRS, 'PRS5':prs5.PRS,'PRS7':prs7.PRS, 'PRS10':prs10.PRS})
prs['SampleID'] = prs1.IID
merged = pd.merge(pheno, prs, on="SampleID").drop(columns=['Duplicate_SUBJID', 'flag_age_covariate'])
merged=merged.rename(columns={"AD_status_final":"Diagnosis", "age_covariate":"Age"})
merged = merged.fillna(-100)
#merged.to_csv('/gpfs/commons/home/tlin/output/prs/bellenguez/fixed_0224/bellenguez_fixed_0224.tsv', sep = '\t', index = False)
merged.to_csv( save + 'kunkle_polypred.tsv', sep = '\t', index = False)
