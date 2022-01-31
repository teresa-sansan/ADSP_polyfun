## phenotype file requires (in the order of) FID IID PHENO
## take SUBJID (#1) SampleID (#2) and AD_status_final(AD diagnosis, #8)
cat /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/phenotype_data_10_28_2021/all_phenotypes_unique_ancestry_all_all.tsv | awk '{OFS = "\t"}{print $1,$2,$8}'> /gpfs/commons/home/tlin/data/ADSP_pheno_merge_beta_new.tsv
