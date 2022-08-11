import pandas as pd
import os

file_name='/gpfs/commons/home/tlin/output/wightman/fixed_0224/susie/finemap_fixed_assertion_susie_iter/polypred/new_beta_max_snp_10_polypred.tsv.prs'
save_name='/gpfs/commons/home/tlin/output/prs/wightman/susie_new_beta_max_snp_10_polypred.tsv'

prs= pd.read_csv(file_name, sep = '\t', names = ["FID","IID","PRS"])
pheno = pd.read_csv("/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/phenotype_data_10_28_2021/all_phenotypes_unique_ancestry_subset.tsv", sep='\t')
pheno_merge = pd.merge(pheno,prs[["IID","PRS"]], left_on ="SampleID", right_on='IID')

pheno_merge.to_csv(save_name,index = False, sep='\t')
