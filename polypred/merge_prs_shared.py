import pandas as pd
from functools import reduce  

path='/gpfs/commons/home/tlin/output/prs/new_anno_0318_24/'
pheno=pd.read_csv('/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/phenotype_file/release_36K/pheno_ADSP_IBD.tsv', sep='\t')

for x in ['bellenguez_adsp_reference']:
    path='/gpfs/commons/home/tlin/output/prs/new_anno_0318_24/'
    name = 'remove_index0'
    path = path + x + '/' + name + '/'
    #without pip thres
    prs_susie = pd.read_csv(path+"susie.prs", sep = ' ', names = ["IID","PRS"])
    prs_bl = pd.read_csv(path+"baseline.prs", sep = ' ', names = ["IID","PRS"])
    prs_bl_omics = pd.read_csv(path+"omics.prs", sep = ' ', names = ["IID","PRS"]) 
    prs_bl_omics_dl = pd.read_csv(path+"omics_dl.prs" , sep = ' ', names = ["IID","PRS"]) 

    prs = [prs_susie, prs_bl, prs_bl_omics, prs_bl_omics_dl]        
    prs_merge = reduce(lambda left, right:pd.merge(left,right,on=["IID"]),prs)
    prs_merge = prs_merge.set_axis( ["IID","PRS_susie","PRS_bl","PRS_bl_omics","PRS_bl_omics_dl"], axis='columns')
    all_merge = pd.merge(pheno,prs_merge, left_on ="SampleID", right_on='IID')
    all_merge.fillna('-1', inplace=True) # Race and Ethniciity
    save= path  + 'prs_pip_%s.tsv'%thres
    all_merge.to_csv(save,index = False, sep='\t')
    print("save prs to " + save)
