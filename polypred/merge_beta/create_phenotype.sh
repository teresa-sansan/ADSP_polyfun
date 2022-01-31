cat all_phenotypes_unique_ancestry_all_all.tsv | awk '{OFS = "\t"}{print $2 ,$1,$8}'> ADSP_pheno_merge_beta.tsv
