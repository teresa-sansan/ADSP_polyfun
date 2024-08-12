import pandas as pd
path="/gpfs/commons/home/tlin/output/sbayesRC/bellenguez_whole_genome/bellenguez_tune"
file_name = '.prs'
#file_name='check_result_bl_19.tsv.prs'

prs= pd.read_csv(path+file_name, sep = ' ', names = ["IID","PRS"])

pheno=pd.read_csv('/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/phenotype_file/release_36K/pheno_ADSP_IBD.tsv', sep='\t')
pheno_merge = pd.merge(pheno,prs[["IID","PRS"]], left_on ="SampleID", right_on='IID')

pheno_merge.to_csv(path+'_pheno.tsv',index = False, sep='\t')
print("done")