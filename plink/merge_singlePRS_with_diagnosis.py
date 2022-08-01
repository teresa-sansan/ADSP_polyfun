import pandas as pd
import os

#pre_file_repo='/gpfs/commons/home/tlin/output/kunkle/kunkle_qc/'
#file_name='/gpfs/commons/home/tlin/output/sbayesR/fixed_0224/sbayesR.prs'
#file_name='/gpfs/commons/home/tlin/output/prs/bellenguez/updateRSID/interested_SNP/maxPIP/updateRSID_interested_SNP_maxPIP.prs'
#save_name='/gpfs/commons/home/tlin/output/prs/bellenguez/updateRSID/interested_SNP/maxPIP/merged_updateRSID_qc_interested_SNP.tsv'

#file_name='/gpfs/commons/home/tlin/output/cT/chr_sep_plink/wightman/new_beta/new_beta_max_snp_10_clump_pT.prs'
#save_name='/gpfs/commons/home/tlin/output/prs/wightman/new_beta_max_snp_10_clump_pT.tsv'

file_name='/gpfs/commons/home/tlin/output/cT/genomewide_plink/kunkle/APOE_qc/kunkle_APOE_qc_chr19.profile'
save_name='/gpfs/commons/home/tlin/output/prs/new_plink/kunkle/APOE_SNP_qc.tsv'

os.system("echo Trimming duplicated spaces...")
no_dup=file_name+"_no_dup_space.prs"
script="cat "+ file_name + "|tr -s ' ' > " + no_dup 
os.system(script)

os.system("echo Merging PRS with phenotype file....")   
#prs= pd.read_csv(no_dup, sep = ' ', names = ["IID","PRS"])
prs= pd.read_csv(no_dup, sep = ' ')
prs= prs.rename(columns={"SCORE": "PRS"})
pheno = pd.read_csv("/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/phenotype_data_10_28_2021/all_phenotypes_unique_ancestry_subset.tsv", sep='\t')
pheno_merge = pd.merge(pheno,prs[["IID","PRS"]], left_on ="SampleID", right_on='IID')

pheno_merge.to_csv(save_name,index = False, sep='\t')
os.system("echo Done!")   
