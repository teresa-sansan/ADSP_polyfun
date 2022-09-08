import pandas as pd

#path = '/gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap/polypred/'
path = '/gpfs/commons/home/tlin/output/bellenguez/bellenguez_all_2/polypred/PLINKupdate_max'
#path = '/gpfs/commons/home/tlin/output/prs/bellenguez_susie/susie_prs_diagnosis_0219.2021_max_snp_'
#path = '/gpfs/commons/home/tlin/output/prs/bellenguez_all_2/with_PC/UPDATEprs_diagnosis_0219.2021_max_snp_'
#path = '/gpfs/commons/home/tlin/output/bellenguez/bellenguez_fixed_0224/finemap/polypred/'
#path = '/gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap/polypred/'
#path = '/gpfs/commons/home/tlin/output/kunkle/kunkle_fixed_0224/finemap/polypred/'
#path='/gpfs/commons/home/tlin/output/wightman/fixed_0224/finemap/polypred/'
#path='/gpfs/commons/home/tlin/output/kunkle/kunkle_fixed_0224/susie_finemap/polypred/'

#path = '/gpfs/commons/home/tlin/output/kunkle/kunkle_fixed_0224_annotations/new_susie/polypred/'

#path='/gpfs/commons/home/tlin/output/kunkle/kunkle_fixed_0224/finemap/polypred_new_plink/'

#path='/gpfs/commons/home/tlin/output/kunkle/kunkle_fixed_0224_annotations/bl/polypred/'

#save = '/gpfs/commons/home/tlin/output/prs/polypred/wightman/fixed_0224.prs.tsv'
#save = '/gpfs/commons/home/tlin/output/prs/polypred/kunkle/susie.prs.tsv'
#name= '_kunkle_susie_polypred.tsv.prs'

#save = '/gpfs/commons/home/tlin/output/prs/polypred/kunkle/new_plink_polypred.tsv'
#save = '/gpfs/commons/home/tlin/output/prs/polypred/kunkle/new_plink_susie.tsv'
#name='_polypred.tsv.prs'


#save = '/gpfs/commons/home/tlin/output/prs/polypred/bellenguez/old_plink_polypred.tsv'
#name='.prs'
 
#save='/gpfs/commons/home/tlin/output/prs/polypred/kunkle/new_plink_bl_polypred.tsv'
#name='_polypred.tsv.prs'

path='/gpfs/commons/home/tlin/output/wightman/fixed_0224_annotations/polypred/susie_max_snp_'
name='_polypred.tsv.prs'
save='/gpfs/commons/home/tlin/output/prs/polypred/wightman/susie'

#prs1 = pd.read_csv(path+'1' + name , sep = '\t')
#prs3 = pd.read_csv(path+'max_snp_3' + name, sep = '\t') 
prs5 = pd.read_csv(path+'5' + name, sep = '\t')
#prs7 = pd.read_csv(path+'max_snp_7' + name, sep = '\t')
prs10 = pd.read_csv(path+'10' + name, sep = '\t')


pheno = pd.read_csv('/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/phenotype_data_10_28_2021/all_phenotypes_unique_ancestry_subset.tsv', sep = '\t')
#prs = pd.DataFrame({'PRS1':prs1.PRS,'PRS3':prs3.PRS, 'PRS5':prs5.PRS,'PRS7':prs7.PRS, 'PRS10':prs10.PRS})
prs = pd.DataFrame({ 'PRS5':prs5.PRS, 'PRS10':prs10.PRS})
#prs = pd.DataFrame({'PRS1':prs1.PRS, 'PRS5':prs5.PRS, 'PRS10':prs10.PRS})
prs['SampleID'] = prs5.IID
merged = pd.merge(pheno, prs, on="SampleID").drop(columns=['Duplicate_SUBJID', 'flag_age_covariate'])
merged=merged.rename(columns={"AD_status_final":"Diagnosis", "age_covariate":"Age"})
merged = merged.fillna(-100)
merged.to_csv( save , sep = '\t', index = False)

print("Finished! Save PRS file to %s"%(save))
