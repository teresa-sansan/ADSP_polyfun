import pandas as pd


name= '_polypred.tsv.prs'
#path='/gpfs/commons/home/tlin/output/wightman/new_anno_0203/'
#path='/gpfs/commons/home/tlin/output/bellenguez/new_anno_0824/'
#path='/gpfs/commons/home/tlin/output/wightman/new_anno_0824_no_partitions/'
path='/gpfs/commons/home/tlin/output/wightman/new_anno_0203/update_all+enformer/finemap/polypred/36k_ibd_max_snp_5_polypred.tsv.prs'
sumstat='wightman'

#anno_list = ["enformer", "no_ml","all_except_enformer", "update_all+enformer", "glasslab"]
#anno_list = ["old_ml"]
#anno_list = ["bl","no_ml","all_enformer","all_anno"]
anno_list = ["update_all+enformer"]
for anno in anno_list:
    print("start running", anno )
    if True:
        # prs1 = pd.read_csv(path+anno+'/finemap/polypred/max_snp_1' + name , sep = '\t')
        # #prs2 = pd.read_csv(path+anno+'/finemap/polypred/max_snp_2' + name, sep = '\t') 
        # prs5 = pd.read_csv(path+anno+'/finemap/polypred/max_snp_5' + name, sep = '\t')
        # #prs7 = pd.read_csv(path+anno+'/finemap/polypred/max_snp_7' + name, sep = '\t')
        # prs10 = pd.read_csv(path+anno+'/finemap/polypred/max_snp_10' + name, sep = '\t')

        #prs1 = pd.read_csv(path+anno+'/finemap/polypred/max_snp_1' + name , sep = '\t')
        #prs3 = pd.read_csv(path+anno+'/finemap/polypred/max_snp_3' + name, sep = '\t')
        #prs5=pd.read_csv(path, sep = '\t')
        prs5 = pd.read_csv(path+anno+'/finemap/polypred/max_snp_5' + name, sep = '\t')
        
        #prs10 = pd.read_csv(path+'/finemap/polypred/max_snp_10' + name, sep = '\t')

        #prs = pd.DataFrame({'PRS1':prs1.PRS,'PRS3':prs3.PRS, 'PRS5':prs5.PRS,'PRS7':prs7.PRS, 'PRS10':prs10.PRS})
        #prs = pd.DataFrame({'PRS1':prs1.PRS, 'PRS5':prs5.PRS,'PRS10':prs10.PRS})   
        prs = pd.DataFrame({ 'PRS5':prs5.PRS}) 
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

    pheno = pd.read_csv('/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/phenotype_file/release_36K/pheno_ADSP_IBD.tsv', sep = '\t')

    merged = pd.merge(pheno, prs, on="SampleID")
    #merged = pd.merge(pheno, prs, on="SampleID").drop(columns=['Duplicate_SUBJID', 'flag_age_covariate'])
    merged = merged.rename(columns={"AD_status_final":"Diagnosis", "age_covariate":"Age"})
    merged = merged.fillna(-100)
    save_name_full = '/gpfs/commons/home/tlin/output/wightman/new_anno_0203/' + sumstat +'.'+anno + '.ibd36k.tsv'
    merged.to_csv(save_name_full, sep = '\t', index = False)
    print("Finished! Save PRS file to %s"%(save_name_full))

#pheno = pd.read_csv('/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/phenotype_data_10_28_2021/all_phenotypes_unique_ancestry_subset.tsv', sep = '\t')