library(ggplot2)
library(lattice)
library(dplyr)
library(stringr)
library(tidyr)
library(ggstance)
library(jtools)
library(tibble)
library(readr)
library(cowplot)
library(pROC)
library(modEvA)  # psuedo rsquare
library(boot)

pre_process <- function(df){
  df<- read.csv(df,sep = '\t', header=T,fill = T)
  if("age_covariate" %in% colnames(df))  ## because there were two files with different column names for phenotype merging.
  {
    colnames(df)[which(names(df) == 'age_covariate')] <- 'Age'
    colnames(df)[which(names(df) == 'AD_status_final')] <- 'Diagnosis'
  }
  df$Age <- as.numeric(as.character(df$Age))
  print(paste("original=",  dim(df)[1], "rows"))
  df <- df %>%
    filter(Diagnosis != -1 & Age >= 65)
  print(paste("filtered=",  dim(df)[1], "rows"))
  return(df)
}

### Load PRS result from different summary statistics and methods. 
### others ----
## new plink
kunkle_adsp <- pre_process('/gpfs/commons/home/tlin/output/prs/new_plink/kunkle/kunkle_ADSP.tsv')
bellenguez_adsp <- pre_process('/gpfs/commons/home/tlin/output/prs/new_plink/bellenguez/bellenguez_ADSP.tsv')
wightman_adsp <- pre_process('/gpfs/commons/home/tlin/output/prs/new_plink/wightman/wightman_ADSP.tsv')


kunkle_adsp_qc <- pre_process('/gpfs/commons/home/tlin/output/prs/new_plink/kunkle/kunkle_ADSP_qc_all.tsv')
bellenguez_adsp_qc <- pre_process('/gpfs/commons/home/tlin/output/prs/new_plink/bellenguez/bellenguez_ADSP_qc_all.tsv')
wightman_adsp_qc <- pre_process('/gpfs/commons/home/tlin/output/prs/new_plink/wightman/wightman_ADSP_qc_all.tsv')


## ADSP phenotype composition ----

pheno_subset <- read.delim("/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/phenotype_data_10_28_2021/all_phenotypes_unique_ancestry_subset.tsv")
pheno_subset %>% count("AD_status_final")
summary(pheno_subset$age_covariate)


### Analysis of QC target file (ADSP plink file)

filtered <- subset(wightman_adsp, !IID %in%  wightman_adsp_qc$IID)
summary(filtered[c('cohort', 'Diagnosis','final_population')])
summary(wightman_adsp[c('cohort', 'Diagnosis','final_population')])

par(mfrow=c(2,3))
draw_pie_chart <- function(df, col, name=' '){
  if (col == 'Diagnosis' || col == 'AD_status_final'){
    df[[col]]<-factor(df[[col]])
    
  }
  pie_labels <- paste(names(summary(df[[col]])),'\n', "(", round(summary(df[[col]])/sum(summary(df[[col]])) * 100,2), "%)")
  pie(summary(df[[col]]), main= name, labels = pie_labels,cex=1,radius=, cex.main=1.5)
}


for (col in c('cohort', 'final_population')){
  draw_pie_chart(wightman_adsp, col, name= 'Before QC')
  draw_pie_chart(filtered, col, name= "Individual been filtered out")
  draw_pie_chart(wightman_adsp_qc, col, name= "After QC")
}

### Check the if the individual being filtered out in genotyping file has duplicates in phenotyping file


remove = read.csv("/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/plink/ADSP_qc_all/ADSP_qc_all.irem", sep = '\t',header = FALSE)
phenotype = read.csv("/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/phenotype_data_10_28_2021/all_phenotypes_unique_ancestry_subset.tsv", sep = '\t')
colnames(remove) <- c('FID','IID')

phenotype_should_remove =phenotype[phenotype$SampleID %in% remove$IID,]
par(mfrow=c(1,3))
draw_pie_chart(phenotype_should_remove, 'Duplicate_SUBJID', name= 'Duplicated SUBJID')
draw_pie_chart(phenotype_should_remove, 'cohort', name= 'Cohort')
draw_pie_chart(phenotype_should_remove, 'AD_status_final', name= 'Diagnosis')
write.table(phenotype_should_remove, file='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/plink/ADSP_qc_all/phenotype_should_remove.tsv', quote = F, sep = '\t')


