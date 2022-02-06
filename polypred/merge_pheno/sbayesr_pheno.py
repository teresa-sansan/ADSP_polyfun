import pandas as pd
prs=pd.read_csv('/gpfs/commons/home/tlin/output/sbayesR/sbayesR.prs', sep=' ', names=["SampleID", 'PRS'])
pheno = pd.read_csv('/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/phenotype_data_10_28_2021/all_phenotypes_unique_ancestry_subset.tsv', sep = '\t')

prs.columns=["SampleID",'PRS']
merged = pheno.merge(prs,on='SampleID')  
merged.to_csv('/gpfs/commons/home/tlin/output/prs/sbayesR.tsv', sep = '\t', index = False)
