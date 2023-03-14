import pandas as pd

#save = '/gpfs/commons/home/tlin/output/prs/polypred/wightman/check_1003_susie.prs.tsv'
#save = '/gpfs/commons/home/tlin/output/prs/polypred/kunkle/susie.prs.tsv'
#save='/gpfs/commons/home/tlin/output/prs/polypred/jansen/new_plink_polypred.tsv '
#save = '/gpfs/commons/home/tlin/output/prs/polypred/bellenguez/susie_polypred.tsv'
#save='/gpfs/commons/home/tlin/output/prs/polypred/bellenguez/fixed_0224_bl_polypred.prs'
save='/gpfs/commons/home/tlin/output/prs/new_plink/wightman/enformer/update_all+enformer.prs'

name= '_polypred.tsv.prs'

#name='.prs'
 
#save='/gpfs/commons/home/tlin/output/prs/polypred/kunkle/new_plink_bl_polypred.tsv'
#name='_polypred.tsv.prs'

# path='/gpfs/commons/home/tlin/output/wightman/fixed_0224_annotations/polypred/susie_max_snp_'
# name='_polypred.tsv.prs'
# save='/gpfs/commons/home/tlin/output/prs/polypred/wightman/susie'

path='/gpfs/commons/home/tlin/output/wightman/wightman_check_1003/bl/finemap/polypred/'
save_name='wightman/check_1003_bl.prs'

if True:
    #prs1 = pd.read_csv(path+'max_snp_1' + name , sep = '\t')
    #prs3 = pd.read_csv(path+'max_snp_3' + name, sep = '\t') 
    prs5 = pd.read_csv(path+'max_snp_5' + name, sep = '\t')
    #prs7 = pd.read_csv(path+'max_snp_7' + name, sep = '\t')
    prs10 = pd.read_csv(path+'max_snp_10' + name, sep = '\t')

    #prs = pd.DataFrame({'PRS1':prs1.PRS,'PRS3':prs3.PRS, 'PRS5':prs5.PRS,'PRS7':prs7.PRS, 'PRS10':prs10.PRS})
   # prs = pd.DataFrame({'PRS1':prs1.PRS, 'PRS5':prs5.PRS, 'PRS10':prs10.PRS})   
    prs = pd.DataFrame({'PRS5':prs5.PRS, 'PRS10':prs10.PRS})   
    prs['SampleID'] = prs5.IID

#path='/gpfs/commons/home/tlin/output/wightman/wightman_check_1003/'
#save_name='wightman/check_1003_fixed_convergence'
#path='/gpfs/commons/home/tlin/output/jansen/'
#save_name='jansen/fixed_convergence'

if False:
    all_anno=pd.read_csv(path+'finemap/polypred/fixed_convergence_max_snp_10_polypred.tsv.prs', sep = '\t')
    #bl = pd.read_csv(path+'bl/finemap/polypred/fixed_max_snp_10_polypred.tsv.prs', sep = '\t')
    susie=pd.read_csv(path+'susie/polypred/fixed_convergence_max_snp_10_polypred.tsv.prs', sep = '\t')
    #prs = pd.DataFrame({'all_anno':all_anno.PRS, 'bl':bl.PRS, 'susie':susie.PRS})   
    prs = pd.DataFrame({'all_anno':all_anno.PRS, 'susie':susie.PRS})   
    prs['SampleID'] = all_anno.IID

pheno = pd.read_csv('/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/phenotype_data_10_28_2021/all_phenotypes_unique_ancestry_subset.tsv', sep = '\t')

merged = pd.merge(pheno, prs, on="SampleID").drop(columns=['Duplicate_SUBJID', 'flag_age_covariate'])
merged = merged.rename(columns={"AD_status_final":"Diagnosis", "age_covariate":"Age"})
merged = merged.fillna(-100)
save_name_full = '/gpfs/commons/home/tlin/output/prs/polypred/' + save_name +'.tsv'
merged.to_csv(save, sep = '\t', index = False)
print("Finished! Save PRS file to %s"%(save_name_full))
