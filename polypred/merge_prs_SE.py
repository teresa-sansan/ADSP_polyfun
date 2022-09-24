import pandas as pd


all_polypred_path = ['/gpfs/commons/home/tlin/output/kunkle/kunkle_fixed_0224/finemap/polypred_new_plink/',
                      '/gpfs/commons/home/tlin/output/bellenguez/old/bellenguez_fixed_0224/finemap/polypred_new_plink/',
                      '/gpfs/commons/home/tlin/output/wightman/wightman_fixed_0224/finemap/polypred_new_plink/',
                     '/gpfs/commons/home/tlin/output/jansen/finemap/polypred/']
                     
pheno = pd.read_csv('/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/phenotype_data_10_28_2021/all_phenotypes_unique_ancestry_subset.tsv', sep = '\t')

name='_polypred.tsv.prs_Mean_SE'

def merge_pheno(path, save_name):
    #merged=merged.rename(columns={"AD_status_final":"Diagnosis", "age_covariate":"Age"})):
    prs1 = pd.read_csv(path+'max_snp_1' + name , sep = '\t')
    #prs3 = pd.read_csv(path+'max_snp_3' + name, sep = '\t') 
    prs5 = pd.read_csv(path+'max_snp_5' + name, sep = '\t')
    #prs7 = pd.read_csv(path+'max_snp_7' + name, sep = '\t')
    prs10 = pd.read_csv(path+'max_snp_10' + name, sep = '\t')

    prs = pd.DataFrame({'PRS1_se':prs1.se, 'PRS1_mean':prs1["mean"], 'PRS5_se':prs5.se, 'PRS5_mean':prs5["mean"], 'PRS10_se':prs10.se,'PRS10_mean':prs10["mean"]})
    prs['SampleID'] = prs5.IID
    merged = pd.merge(pheno, prs, on="SampleID").drop(columns=['Duplicate_SUBJID', 'flag_age_covariate'])
    merged=merged.rename(columns={"AD_status_final":"Diagnosis", "age_covariate":"Age"})
    merged.to_csv(save_name, sep = '\t', index = False)
    print("Finished! Save PRS file to %s"%(save_name))


kunkle = merge_pheno(all_polypred_path[0],'/gpfs/commons/home/tlin/output/prs/sumstat_jk/kunkle_new_plink_jk.tsv')
bellenguez = merge_pheno(all_polypred_path[1],'/gpfs/commons/home/tlin/output/prs/sumstat_jk/bellenguez_new_plink_jk.tsv')
wightman = merge_pheno(all_polypred_path[2],'/gpfs/commons/home/tlin/output/prs/sumstat_jk/wightman_new_plink_jk.tsv')