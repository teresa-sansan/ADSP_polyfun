import pandas as pd
import os

#file_name='/gpfs/commons/home/tlin/output/wightman/fixed_0224/susie/finemap_fixed_assertion_susie_iter/polypred/new_beta_max_snp_10_polypred.tsv.prs'
#file_name='/gpfs/commons/home/tlin/output/bellenguez/bellenguez_fixed_0224_annotations/susie/max_snp_10/polypred/susie_max_snp_10_polypred.tsv.prs'
#save_name='/gpfs/commons/home/tlin/output/prs/polypred/bellenguez/susie_max_snp_10_polypred.tsv'

#file_name='/gpfs/commons/home/tlin/output/jansen/finemap/polypred/max_snp_10_polypred.tsv.prs'
#save_name='/gpfs/commons/home/tlin/output/prs/polypred/jansen/max_snp_10_polypred.tsv.prs'

#file_name='/gpfs/commons/home/tlin/output/bellenguez/old/bellenguez_fixed_0224/finemap/polypred_new_plink/not_fixed_max_snp_10_polypred.tsv.prs'
#save_name='/gpfs/commons/home/tlin/output/bellenguez/old/bellenguez_fixed_0224/finemap/polypred_new_plink/max_snp_10_no_rescue'

# file_name='/gpfs/commons/home/tlin/output/wightman/wightman_check_1003/susie/finemap/polypred/fixed_max_snp_10_polypred.tsv.prs'
# save_name='/gpfs/commons/home/tlin/output/prs/polypred/wightman/check_1003_susie_max_snp_10_rescue'
file_name='/gpfs/commons/home/tlin/output/wightman/new_anno_0203/update_all+enformer/finemap/polypred/w_prscs/polypred.predictions.prs'
save_name='/gpfs/commons/home/tlin/output/wightman/new_anno_0203/update_all+enformer/finemap/polypred/w_prscs/polypred.prs.wpheno'

prs= pd.read_csv(file_name, sep = '\t', names = ["IID","FID","PRS"])
pheno = pd.read_csv("/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/phenotype_data_10_28_2021/all_phenotypes_unique_ancestry_subset.tsv", sep='\t')
pheno_merge = pd.merge(pheno,prs, left_on ="SampleID", right_on='IID')
pheno_merge.to_csv(save_name+'.tsv',index = False, sep='\t')

# file_name=['/gpfs/commons/home/tlin/output/kunkle/kunkle_fixed_0224/finemap/polypred/10.prs_Mean_SE',
#           '/gpfs/commons/home/tlin/output/bellenguez/bellenguez_fixed_0224/finemap/max_snp_10/try_rescue_not_converge/polypred/polypred_genomewide.tsv.prs_Mean_SE',
#           '/gpfs/commons/home/tlin/output/wightman/wightman_fixed_0224/finemap/polypred/max_snp_10_new_beta_polypred.tsv.prs_Mean_SE',
#           '/gpfs/commons/home/tlin/output/jansen/finemap/polypred/max_snp_10_polypred.tsv.prs_Mean_SE']

# save_name='/gpfs/commons/home/tlin/output/prs/sumstat_jk/'
# sumstat=['kunkle.tsv', 'bellenguez.tsv','wightman.tsv','jansen.tsv']
# for i in range(0,4):
#     prs= pd.read_csv(file_name[i], sep = '\t', names = ["IID","FID","prs_mean","prs_sd"])
#     pheno = pd.read_csv("/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/phenotype_data_10_28_2021/all_phenotypes_unique_ancestry_subset.tsv", sep='\t')
#     pheno_merge = pd.merge(pheno,prs, left_on ="SampleID", right_on='IID')
#     pheno_merge.to_csv(save_name+sumstat[i],index = False, sep='\t')


    
    
# file_name='/gpfs/commons/home/tlin/output/jansen/finemap/polypred/max_snp_10_polypred.tsv.prs'
# save_name='/gpfs/commons/home/tlin/output/prs/polypred/jansen/max_snp_10_polypred.tsv.prs'
# prs= pd.read_csv(file_name, sep = '\t', names = ["FID","IID","PRS"])
# pheno = pd.read_csv("/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/phenotype_data_10_28_2021/all_phenotypes_unique_ancestry_subset.tsv", sep='\t')
# pheno_merge = pd.merge(pheno,prs[["IID","PRS"]], left_on ="SampleID", right_on='IID')

# pheno_merge.to_csv(save_name,index = False, sep='\t')