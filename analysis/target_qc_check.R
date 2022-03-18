if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("snpStats")
install.packages("snpStats")
library(snpStats)

chr22_qc = read.plink('/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/plink_biallelic/qc/ADSP_qc_chr22')
chr22 = read.plink('/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/plink_biallelic/ADSP_chr22')

length(chr22$fam$member)     ##16906
length(chr22_qc$fam$member)  ##10157


all_pheno <- read.table('/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/phenotype_data_10_28_2021/all_phenotypes_unique_ancestry_subset.tsv', header =T, fill = T)
not_seen = chr22$fam$member[!(chr22$fam$member %in% chr22_qc$fam$member)] ##6749

not_seen_pheno <-pheno_all %>%                  ##5132
  filter(pheno_all$SampleID %in% not_seen) %>%  ##4186
  filter(!age_covariate < 65)


summarystat <- load('/gpfs/commons/home/tlin/susie/SummaryConsistency1k.RData')

## double check
chr1_qc = read.plink('/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/plink_biallelic/qc/ADSP_qc_chr1')
chr1 = read.plink('/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/plink_biallelic/filt/1_filt') ## there are duplicates so I can't read the file 

