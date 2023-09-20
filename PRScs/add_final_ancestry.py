import pandas as pd

prs = pd.read_csv('prs_36k.tsv', sep = '\t')
pheno = pd.read_csv('/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/phenotype_file/release_36K/pheno_LOAD_1000k.tsv', sep = '\t')
pheno_sub = pheno.iloc[:, list([1,10]) + list(range(12, 32))]

prs = pd.merge(prs,pheno_sub, on='SampleID')
prs = prs.rename(columns={'predicted_ancestry': 'final_population'})
prs.to_csv('prscs_new_ancestry.tsv', sep = '\t', index = False) 