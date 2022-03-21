import pandas as pd
import os


#pre_file_repo='/gpfs/commons/home/tlin/output/kunkle/kunkle_qc/'
#file_name='/gpfs/commons/home/tlin/output/sbayesR/fixed_0224/sbayesR.prs'
file_name='/gpfs/commons/home/tlin/output/prs/bellenguez/updateRSID/interested_SNP/updateRSID_qc_interested_SNP.prs'
save_name='/gpfs/commons/home/tlin/output/prs/bellenguez/updateRSID/interested_SNP/merged_updateRSID_qc_interested_SNP.tsv'


#os.system("cd /gpfs/commons/home/tlin/output/kunkle/kunkle_qc/")
#os.system("pwd")
#os.system("echo trimming duplicated spaces!")
#os.system("cat updated_0224_2SNP_qc.prs.profile|tr -s ' ' > fixed_0224_no_dup_space.prs")


prs= pd.read_csv(file_name, sep = ' ', names = ["IID","PRS"])

pheno = pd.read_csv("/gpfs/commons/home/tlin/output/prs/bellenguez_all_2/with_PC/UPDATEprs_diagnosis_0219.2021_max_snp_10_subset.tsv", sep='\t')
pheno_merge = pd.merge(pheno.drop(columns = "PRS"),prs[["IID","PRS"]], left_on ="SampleID", right_on='IID')

pheno_merge.to_csv(save_name,index = False, sep='\t')
