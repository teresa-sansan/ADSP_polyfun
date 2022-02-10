import pandas as pd
prs = pd.read_csv('/gpfs/commons/home/tlin/output/prs/PRSice/PRSice_pheno.best', sep = ' ')

#prs=pd.read_csv('/gpfs/commons/home/tlin/output/sbayesR/sbayesR.prs', sep=' ', names=["SampleID", 'PRS'])
#prs=pd.read_csv('/gpfs/commons/home/tlin/output/prs/bellenguez_susie/bellenguez_susie_prs.tsv', sep='\t')
#prs = pd.read_csv('/gpfs/commons/home/tlin/output/prs/bellenguez_updateRSID_CS.tsv', sep = '\t')

pheno = pd.read_csv('/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/phenotype_data_10_28_2021/all_phenotypes_unique_ancestry_subset.tsv', sep = '\t')

#prs.columns=["SampleID",'PRS']
#merged = pheno.merge(prs,on='SampleID')  
#merged = pheno.merge(prs.loc[:,["SampleID","PRS1","PRS3","PRS5","PRS7","PRS10"]],left_on='SampleID', right_on='IID')  
#merged = pheno.merge(prs.loc[:,["IID","PRS1","PRS3","PRS5","PRS7","PRS10"]],left_on='SampleID', right_on='IID')  

merged = pheno.merge(prs, left_on='SampleID', right_on='IID')
merged.to_csv('/gpfs/commons/home/tlin/output/prs/PRSice_pheno.tsv', sep = '\t', index = False)

#merged.to_csv('/gpfs/commons/home/tlin/output/prs/bellenguez_updateRSID_prs_PC.tsv', sep = '\t', index = False)
#merged.to_csv('/gpfs/commons/home/tlin/output/prs/bellenguez_susie/bellenguez_susie_prs_PC.tsv', sep = '\t', index = False)
#merged.to_csv('/gpfs/commons/home/tlin/output/prs/sbayesR.tsv', sep = '\t', index = False)
