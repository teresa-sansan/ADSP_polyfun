import pandas as pd
from functools import reduce  

#path='/gpfs/commons/home/tlin/output/prs/pT_36k_ibd/'
path='/gpfs/commons/home/tlin/output/prs/new_anno_0318_24/'
pheno=pd.read_csv('/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/phenotype_file/release_36K/pheno_ADSP_IBD.tsv', sep='\t')
#/gpfs/commons/home/tlin/output/prs/new_anno_0318_24/bellenguez_adsp_reference/allele_flip//omics_dl.prs

for x in ['bellenguez_adsp_reference']:
    path='/gpfs/commons/home/tlin/output/prs/new_anno_0318_24/'
    # path= path + x + '/thres/plink_output/'
    path= path + x + '/allele_flip/'
    
    ##without thres
    prs_susie = pd.read_csv(path+"susie.prs", sep = ' ', names = ["IID","PRS"])
    prs_bl = pd.read_csv(path+"baseline.prs", sep = ' ', names = ["IID","PRS"])
    prs_bl_omics = pd.read_csv(path+"omics.prs", sep = ' ', names = ["IID","PRS"]) 
    prs_bl_omics_dl = pd.read_csv(path+"omics_dl.prs" , sep = ' ', names = ["IID","PRS"]) 
    #for thres in [0.1,0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8 ]:

    
        # prs_susie = pd.read_csv(path+"sum_susie_pip_%s.prs"%thres, sep = ' ', names = ["IID","PRS"])
        # prs_bl = pd.read_csv(path+"sum_baseline_pip_%s.prs"%thres, sep = ' ', names = ["IID","PRS"])
        # prs_bl_omics = pd.read_csv(path+"sum_omics_pip_%s.prs"%thres, sep = ' ', names = ["IID","PRS"]) 
        # prs_bl_omics_dl = pd.read_csv(path+"sum_omics_dl_pip_%s.prs"%thres, sep = ' ', names = ["IID","PRS"]) 

    prs = [prs_susie, prs_bl, prs_bl_omics, prs_bl_omics_dl]        
    prs_merge = reduce(lambda left, right:pd.merge(left,right,on=["IID"]),prs)
    prs_merge = prs_merge.set_axis( ["IID","PRS_susie","PRS_bl","PRS_bl_omics","PRS_bl_omics_dl"], axis='columns')
    all_merge = pd.merge(pheno,prs_merge, left_on ="SampleID", right_on='IID')
    all_merge.fillna('-1', inplace=True) # Race and Ethniciity
    #save= path  + 'prs_pip_%s.tsv'%thres
    save= path  + 'allele_flip.prs.tsv'
    all_merge.to_csv(save,index = False, sep='\t')
    print("save prs to " + save)
