import pandas as pd

#path='/gpfs/commons/home/tlin/output/cT/bellenguez/fixed_0224/qc_on_variant/'
#save_name='bellenguez/fixed_0224/bellenguez_qc_on_variant'

#path='/gpfs/commons/home/tlin/output/cT/kunkle/fixed_0224/qc_on_variant_sumstat/'
#save_name='kunkle/fixed_0224/qc_on_variant_sumstat'

##kunkle_no_apoe
#path='/gpfs/commons/home/tlin/output/cT/kunkle/fixed_0224/subsets/qc_on_variant_sumstat/'
#save_name='kunkle/fixed_0224/remove_APOE_qc_on_variant_sumstat.tsv'
#path='/gpfs/commons/home/tlin/output/cT/genomewide_plink/kunkle/ADSP_no_apoe/'
#save_name='/gpfs/commons/home/tlin/output/prs/new_plink/kunkle/kunkle_no_apoe_qc'


#path='/gpfs/commons/home/tlin/output/cT/wightman/before_qc/'
#save_name='kunkle/fixed_0224/qc_on_variant_maf01'
#path='/gpfs/commons/home/tlin/output/cT/wightman/qc_on_variant_sumstat/'
#save_name='wightman/qc_on_variant_sumstat'

#path='/gpfs/commons/home/tlin/output/cT/bellenguez/fixed_0224/qc_on_variant_sumstat/'
#save_name='bellenguez/fixed_0224/bellenguez_qc_on_variant_sumstat'

## new beta
#path='/gpfs/commons/home/tlin/output/cT/kunkle/fixed_0224/new_beta/'
#save_name='kunkle/fixed_0224/new_beta_noqc'


sumstat=["bellenguez","kunkle","wightman"]
#plink=["ADSP","ADSP_qc_all","ADSP_qc_variant","ADSP_UKBB","ADSP_UKBB_qc"]
#sumstat=["kunkle"]
plink=["ADSP_no_apoe"]

for x in sumstat:
  #for y in plink:
  #path = '/gpfs/commons/home/tlin/output/cT/genomewide_plink/'+ x + '/' + y + '/'
  path = '/gpfs/commons/home/tlin/output/cT/new_plink/' + x + '/fixed_0224/polyfun_beta/'
  save_name = 'new_plink/' + x + '/' + x + "_" + 'polyfun_beta' 
  prs_e5 = pd.read_csv(path+"pT_e-5.prs", sep = ' ', names = ["IID","PRS"])
  prs_001 = pd.read_csv(path+"pT_0.001.prs", sep = ' ', names = ["IID","PRS"])
  prs_005 = pd.read_csv(path+"pT_0.005.prs", sep = ' ', names = ["IID","PRS"])   
  prs_01 = pd.read_csv(path+"pT_0.01.prs", sep = ' ', names = ["IID","PRS"]) 
  prs_05 = pd.read_csv(path+"pT_0.05.prs", sep = ' ', names = ["IID","PRS"])
  prs_1 = pd.read_csv(path+"pT_0.1.prs", sep = ' ', names = ["IID","PRS"]) 
  prs_5 = pd.read_csv(path+"pT_0.5.prs", sep = ' ', names = ["IID","PRS"]) 

  pheno = pd.read_csv("/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/phenotype_data_10_28_2021/all_phenotypes_unique_ancestry_subset.tsv", sep='\t')
  prs = [prs_e5, prs_001, prs_005, prs_01, prs_05, prs_1, prs_5]  
  from functools import reduce  
  prs_merge = reduce(lambda left, right:pd.merge(left,right,on=["IID"]),prs)

  prs_merge = prs_merge.set_axis( ["IID","PRS_e5","PRS_001","PRS_005","PRS_01","PRS_05","PRS_1","PRS_5"], axis='columns')
  all_merge = pd.merge(pheno,prs_merge, left_on ="SampleID", right_on='IID')
  all_merge.fillna('-1', inplace=True) # Race and Ethniciity
  save="/gpfs/commons/home/tlin/output/prs/"+ save_name+'.tsv'
  all_merge.to_csv(save,index = False, sep='\t')
  print("save prs to " + save)
