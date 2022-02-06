library(ggplot2)
library(lattice)
library(dplyr)
library(stringr)
library(tidyr)
library(ggstance)
library(jtools)
library(tibble)
library(readr)


## first, look at the composition of ADSP phenotype file (in plink format) ----

pT <- pre_process("/gpfs/commons/home/tlin/output/prs/bellenguez_pT_PRS_withPC.tsv")
pheno_all <- read.delim("/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/phenotype_data_10_28_2021/all_phenotypes_unique_ancestry_all_all.tsv")
pheno_subset <- read.delim("/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/phenotype_data_10_28_2021/all_phenotypes_unique_ancestry_subset.tsv")

pheno_all %>% count(AD_status_final)
pheno_subset %>% count(AD_status_final)
pheno_subset %>% count(age_covariate)
pT %>% count(final_population)
pT[pT$final_population=="EUR",]%>% count(Diagnosis)


## load PRS data ----
sbayesR <- read.delim("/gpfs/commons/home/tlin/output/prs/sbayesR.tsv")  
names(sbayesR)[names(sbayesR) == 'age_covariate'] <- "Age"
names(sbayesR)[names(sbayesR) == 'AD_status_final'] <- "Diagnosis"
sbayesR = pre_process(sbayesR, FILE = TRUE)
  ##13265 -> 11721

density_plot <- function(polypred, name){
  sort = polypred[order(polypred$PRS),]$PRS
  plot(density(polypred[polypred$Diagnosis == 1,]$PRS),col = "red", main=name, xlab="PRS") 
  lines(density(polypred[polypred$Diagnosis == 0,]$PRS), col = "blue")
  legend("topright", legend=c("Case", "Control"),
         col=c("red", "blue"), lty=1:1, cex=0.8,
         box.lty=0)
}

density_plot(sbayesR, name = "sbayesR")

par(mfrow=c(2,3))

p_threshold_col <- c("PRS_005","PRS_05","PRS_1","PRS_5")
pT_subset = pT[p_threshold_col]

for(i in c(33,35,37,34,36,38)) {       # for-loop over columns
  plot(density((pT[, i])),col='grey', lty = 2, lwd = 2,main=str_replace(colnames(pT)[i],"PRS_", "p=0."), xlab='PRS')
  #plot(density(pT[pT$Diagnosis == 1,i]), col = 'red',main=str_replace(colnames(pT)[i],"PRS_", "p = 0."), xlab='PRS')
  lines(density(pT[pT$Diagnosis == 1,i]), col = 'red')
  lines(density(pT[pT$Diagnosis == 0,i]), col = "blue")
  legend("topright", legend=c("ALL","Case", "Control"),
         col=c("darkgrey","red", "blue"), lty=c(2,1,1), cex=0.8,
         box.lty=0)
}




