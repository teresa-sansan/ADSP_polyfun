import pandas as pd
from functools import reduce  

#path='/gpfs/commons/home/tlin/output/cT/new_plink_genomewide/wightman/fixed_rsid_1002/ADSP_qc_all/'
#save_name='wightman/fixed_rsid1002/ADSP_qc_all'

# new beta for wightman

#sumstat=["bellenguez","kunkle","wightman"]
#sumstat=['wightman']

#plink=["ADSP","ADSP_qc_all","ADSP_qc_variant","ADSP_UKBB","ADSP_UKBB_qc"]
#plink=["ADSP_no_apoe"]
#path = '/gpfs/commons/home/tlin/output/cT/new_plink_genomewide/kunkle/ADSP_no_apoe/'

#path='/gpfs/commons/home/tlin/output/cT/new_plink_genomewide/bellenguez/new_sep22/ADSP/'
#save_name='bellenguez/new_sep22/ADSP'


#path='/gpfs/commons/home/tlin/output/cT/new_plink_genomewide/jansen/ADSP/'
#save_name='jansen/ADSP'


#sumstat=["jansen"]
#plink=["ADSP_qc_all"]
#sumstat=['kunkle']

#path='/gpfs/commons/home/tlin/output/prs/pT_36k_ibd/'
path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_QC/annotated_hg37_plink_ibd/plink_prs/clump_bellenguez/'
save_name='bellenguez_pT_36k_ibd'



## Polyfun thres
#path='/gpfs/commons/home/tlin/output/prs/new_anno_0318_24/bellenguez_adsp_reference/thres/plink_output/'
path='/gpfs/commons/home/tlin/output/prs/new_anno_0318_24/bellenguez_adsp_reference/remove_index0/'
## susie doesht have p < 0.0006
for anno in ["susie","baseline","omics","omics_dl"]:
    
    prs_1 = pd.read_csv(path+ anno+"_pip_0.1.prs", sep = ' ', names = ["IID","PRS"])
    prs_2 = pd.read_csv(path+ anno+"_pip_0.2.prs", sep = ' ', names = ["IID","PRS"])
    prs_3 = pd.read_csv(path+ anno+"_pip_0.3.prs", sep = ' ', names = ["IID","PRS"])
    prs_4 = pd.read_csv(path+ anno+"_pip_0.4.prs", sep = ' ', names = ["IID","PRS"])
    prs_5 = pd.read_csv(path+ anno+"_pip_0.5.prs", sep = ' ', names = ["IID","PRS"])
    prs_6 = pd.read_csv(path+ anno+"_pip_0.6.prs", sep = ' ', names = ["IID","PRS"])
    prs_7 = pd.read_csv(path+ anno+"_pip_0.7.prs", sep = ' ', names = ["IID","PRS"])
    prs_8 = pd.read_csv(path+ anno+"_pip_0.8.prs", sep = ' ', names = ["IID","PRS"])
    
    pheno=pd.read_csv('/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/phenotype_file/release_36K/pheno_ADSP_IBD.tsv', sep='\t')
    prs = [prs_1,prs_2,prs_3,prs_4,prs_5,prs_6,prs_7,prs_8]  
    prs_merge = reduce(lambda left, right:pd.merge(left,right,on=["IID"]),prs)
    prs_merge = prs_merge.set_axis( ["IID","PRS_1","PRS_2","PRS_3","PRS_4","PRS_5","PRS_6","PRS_7","PRS_8"], axis='columns')

    
    # prs_e4 = pd.read_csv(path+anno+"_e-4prs", sep = ' ', names = ["IID","PRS"])
    # prs_e5 = pd.read_csv(path+anno+"_e-5prs", sep = ' ', names = ["IID","PRS"])
    # #prs_e6 = pd.read_csv(path+anno+"_e-6prs", sep = ' ', names = ["IID","PRS"])
    # prs_001 = pd.read_csv(path+anno+"_0.001prs", sep = ' ', names = ["IID","PRS"])
    # prs_01 = pd.read_csv(path+anno+"_0.01prs", sep = ' ', names = ["IID","PRS"]) 
    # prs_1 = pd.read_csv(path+anno+"_0.1prs", sep = ' ', names = ["IID","PRS"]) 
    # pheno=pd.read_csv('/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/phenotype_file/release_36K/pheno_ADSP_IBD.tsv', sep='\t')
    # prs = [prs_e5, prs_e4, prs_001, prs_01,  prs_1]  
    # prs_merge = reduce(lambda left, right:pd.merge(left,right,on=["IID"]),prs)
    # prs_merge = prs_merge.set_axis( ["IID","PRS_e5","PRS_e4","PRS_001","PRS_01","PRS_1"], axis='columns')
    #prs_merge = prs_merge.set_axis( ["IID","PRS_e5","PRS_001","PRS_005","PRS_01","PRS_05","PRS_1","PRS_5"], axis='columns')
    all_merge = pd.merge(pheno,prs_merge, left_on ="SampleID", right_on='IID')
    all_merge.fillna('-1', inplace=True) # Race and Ethniciity
    all_merge.to_csv('/gpfs/commons/home/tlin/output/prs/new_anno_0318_24/bellenguez_adsp_reference/thres/remove_index0_'+anno+'.tsv',index = False, sep='\t')
    #print(path+anno+'.tsv')

# prs_e4 = pd.read_csv(path+"pT_e-4.prs", sep = ' ', names = ["IID","PRS"])
# prs_e5 = pd.read_csv(path+"pT_e-5.prs", sep = ' ', names = ["IID","PRS"])
# prs_e6 = pd.read_csv(path+"pT_e-6.prs", sep = ' ', names = ["IID","PRS"])
# prs_001 = pd.read_csv(path+"pT_0.001.prs", sep = ' ', names = ["IID","PRS"])
# #prs_005 = pd.read_csv(path+"pT_0.005.prs", sep = ' ', names = ["IID","PRS"])   
# prs_01 = pd.read_csv(path+"pT_0.01.prs", sep = ' ', names = ["IID","PRS"]) 
# #prs_05 = pd.read_csv(path+"pT_0.05.prs", sep = ' ', names = ["IID","PRS"])
# prs_1 = pd.read_csv(path+"pT_0.1.prs", sep = ' ', names = ["IID","PRS"]) 
# #prs_5 = pd.read_csv(path+"pT_0.5.prs", sep = ' ', names = ["IID","PRS"]) 

# #pheno = pd.read_csv("/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/phenotype_data_10_28_2021/all_phenotypes_unique_ancestry_subset.tsv", sep='\t')
# #pheno=pd.read_csv("/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/phenotype_file/release_36K/pheno_LOAD_1000k.tsv", sep = '\t')
# pheno=pd.read_csv('/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/phenotype_file/release_36K/pheno_ADSP_IBD.tsv', sep='\t')
# #prs = [prs_e5, prs_001, prs_005, prs_01, prs_05, prs_1, prs_5]  
# prs = [prs_e6,prs_e5, prs_e4, prs_001, prs_01,  prs_1]  
# prs_merge = reduce(lambda left, right:pd.merge(left,right,on=["IID"]),prs)
# prs_merge = prs_merge.set_axis( ["IID","PRS_e6","PRS_e5","PRS_e4","PRS_001","PRS_01","PRS_1"], axis='columns')
# #prs_merge = prs_merge.set_axis( ["IID","PRS_e5","PRS_001","PRS_005","PRS_01","PRS_05","PRS_1","PRS_5"], axis='columns')
# all_merge = pd.merge(pheno,prs_merge, left_on ="SampleID", right_on='IID')
# all_merge.fillna('-1', inplace=True) # Race and Ethniciity
# save="/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/PRS/36k_ibd/"+ save_name+'.tsv'
# all_merge.to_csv(save,index = False, sep='\t')
# print("save prs to " + save)
