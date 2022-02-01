install pandas as pd
prs=pd.read_csv('/gpfs/commons/home/tlin/output/sbayesR/polypred/polypred.pred.prs', sep='\t')
pheno = pd.read_csv('/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/phenotype_data_10_28_2021/all_phenotypes_unique_ancestry_subset.tsv', sep = '\t')
merged = pheno.merge(prs,left_on='SampleID', right_on='IID').drop(columns = 'FID')  
merged.to_csv('/gpfs/commons/home/tlin/output/prs/sbayesr_polyfun/prs_pheno.tsv', sep = '\t', index = False)
