library(ggplot2)
library(lattice)
library(dplyr)
library(stringr)
library(tidyr)
library(ggstance)
library(jtools)
library(tibble)
library(readr)


## first, look at the composition of ADSP phenotype file (in plinl format)

pT <- read.delim("/gpfs/commons/home/tlin/output/prs/bellenguez_all_2/pT/pT_PRS_withPC.tsv")
pheno_all <- read.delim("/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/phenotype_data_10_28_2021/all_phenotypes_unique_ancestry_all_all.tsv")
pheno_subset <- read.delim("/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/phenotype_data_10_28_2021/all_phenotypes_unique_ancestry_subset.tsv")
pT <- pre_process("/gpfs/commons/home/tlin/output/prs/bellenguez_all_2/pT/pT_PRS_withPC.tsv")
pheno_all %>% count(AD_status_final)
pheno_subset %>% count(AD_status_final)
pheno_subset %>% count(age_covariate)
pT %>% count(final_population)
pT[pT$final_population=="EUR",]%>% count(Diagnosis)



updatePRS1 = pre_process("/gpfs/commons/home/tlin/output/prs/bellenguez_all_2/with_PC/UPDATEprs_diagnosis_0219.2021_max_snp_1_subset.tsv")
updatePRS3 = pre_process("/gpfs/commons/home/tlin/output/prs/bellenguez_all_2/with_PC/UPDATEprs_diagnosis_0219.2021_max_snp_3_subset.tsv")
updatePRS5 = pre_process("/gpfs/commons/home/tlin/output/prs/bellenguez_all_2/with_PC/UPDATEprs_diagnosis_0219.2021_max_snp_5_subset.tsv")
updatePRS7 = pre_process("/gpfs/commons/home/tlin/output/prs/bellenguez_all_2/with_PC/UPDATEprs_diagnosis_0219.2021_max_snp_7_subset.tsv")
updatePRS10 = pre_process("/gpfs/commons/home/tlin/output/prs/bellenguez_all_2/with_PC/UPDATEprs_diagnosis_0219.2021_max_snp_10_subset_new.tsv")

