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

path='/gpfs/commons/home/tlin/output/prs/pT_36k/wightman/'
save_name='wightman_pT'

#for x in sumstat:
  #for y in plink:
  #path = '/gpfs/commons/home/tlin/output/cT/old_plink_chr_sep/'+ x + '/' + y + '/'
  #print("start " + x)
prs_e5 = pd.read_csv(path+"pT_e-5.prs", sep = ' ', names = ["IID","PRS"])
prs_001 = pd.read_csv(path+"pT_0.001.prs", sep = ' ', names = ["IID","PRS"])
prs_005 = pd.read_csv(path+"pT_0.005.prs", sep = ' ', names = ["IID","PRS"])   
prs_01 = pd.read_csv(path+"pT_0.01.prs", sep = ' ', names = ["IID","PRS"]) 
prs_05 = pd.read_csv(path+"pT_0.05.prs", sep = ' ', names = ["IID","PRS"])
prs_1 = pd.read_csv(path+"pT_0.1.prs", sep = ' ', names = ["IID","PRS"]) 
prs_5 = pd.read_csv(path+"pT_0.5.prs", sep = ' ', names = ["IID","PRS"]) 

#pheno = pd.read_csv("/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/phenotype_data_10_28_2021/all_phenotypes_unique_ancestry_subset.tsv", sep='\t')
pheno=pd.read_csv("/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/phenotype_file/release_36K/pheno_LOAD_1000k.tsv", sep = '\t')
prs = [prs_e5, prs_001, prs_005, prs_01, prs_05, prs_1, prs_5]  
prs_merge = reduce(lambda left, right:pd.merge(left,right,on=["IID"]),prs)

prs_merge = prs_merge.set_axis( ["IID","PRS_e5","PRS_001","PRS_005","PRS_01","PRS_05","PRS_1","PRS_5"], axis='columns')
all_merge = pd.merge(pheno,prs_merge, left_on ="SampleID", right_on='IID')
all_merge.fillna('-1', inplace=True) # Race and Ethniciity
save="/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/PRS/36k/"+ save_name+'.tsv'
all_merge.to_csv(save,index = False, sep='\t')
print("save prs to " + save)
