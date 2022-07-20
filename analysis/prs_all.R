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
#library(fmsb)   ##psudo rsquare
library(modEvA)  # psuedo rsquare
library(boot)
library("gridGraphics")

## first, look at the composition of ADSP phenotype file (in plink format) ----

pT <- pre_process("/gpfs/commons/home/tlin/output/prs/bellenguez_pT_PRS_withPC.tsv")
pheno_all <- read.delim("/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/phenotype_data_10_28_2021/all_phenotypes_unique_ancestry_all_all.tsv")
pheno_subset <- read.delim("/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/phenotype_data_10_28_2021/all_phenotypes_unique_ancestry_subset.tsv")

pheno_all %>% count(AD_status_final)
pheno_subset %>% count(AD_status_final)
pheno_subset %>% count(age_covariate)
pT %>% count(final_population)
pT[pT$final_population=="EUR",]%>% count(Diagnosis)

pre_process <- function(df, FILE=FALSE){
  if(FILE==FALSE){
    df<- read.csv(df,sep = '\t', header=T,fill = T)
  }
  if("age_covariate" %in% colnames(df))
  {
    print("change column name first...")
    colnames(df)[which(names(df) == 'age_covariate')] <- 'Age'
    colnames(df)[which(names(df) == 'AD_status_final')] <- 'Diagnosis'
  }
  df$Age <- as.numeric(as.character(df$Age))
  #df$Age <- as.numeric(df$Age)
  print(paste("original=",  dim(df)[1], "rows"))
  
  df <- df %>%
    filter(Diagnosis != -1 & Age >= 65)
  print(paste("filtered=",  dim(df)[1], "rows"))
  return(df)
}


## load PRS data ----

### bellenguez ----- 
bellenguez_updateRSID <- pre_process('/gpfs/commons/home/tlin/output/prs/bellenguez/updateRSID/bellenguez_polypred_prs.tsv')
bellenguez_fixed_0224 <- pre_process("/gpfs/commons/home/tlin/output/prs/bellenguez/fixed_0224/bellenguez_fixed_0224.tsv")
bellenguez_cT <- pre_process('/gpfs/commons/home/tlin/output/prs/bellenguez/bellenguez_pT_PRS_withPC.tsv')
bellenguez_qc <- pre_process('/gpfs/commons/home/tlin/output/prs/bellenguez/updateRSID/bellenguez_qc_pT_PRS.tsv')
bellenguez_qc_on_target <- pre_process('/gpfs/commons/home/tlin/output/prs/bellenguez/updateRSID/qc_on_target.tsv')
bellenguez_qc_on_base <- pre_process('/gpfs/commons/home/tlin/output/prs/bellenguez/updateRSID/qc_on_base.tsv')
bellenguez_cT0224 <- pre_process('/gpfs/commons/home/tlin/output/prs/bellenguez/fixed_0224/bellenguez_no_qc.tsv')
bellenguez_qc0224 <- pre_process('/gpfs/commons/home/tlin/output/prs/bellenguez/fixed_0224/bellenguez_qc.tsv')
bellenguez_qc_on_target0224 <- pre_process('/gpfs/commons/home/tlin/output/prs/bellenguez/fixed_0224/bellenguez_qc_on_target.tsv')
bellenguez_qc_on_base0224 <- pre_process('/gpfs/commons/home/tlin/output/prs/bellenguez/fixed_0224/bellenguez_qc_on_base.tsv')
bellenguez_interested <- pre_process('/gpfs/commons/home/tlin/output/prs/bellenguez/updateRSID/interested_SNP/merged_updateRSID_interested_SNP.tsv') 
bellenguez_qc_interested <- pre_process('/gpfs/commons/home/tlin/output/prs/bellenguez/updateRSID/interested_SNP/merged_updateRSID_qc_interested_SNP.tsv')
bellenguez_interest_max <- pre_process('/gpfs/commons/home/tlin/output/prs/bellenguez/updateRSID/interested_SNP/maxPIP/merged_updateRSID_qc_interested_SNP.tsv')
bellenguez_interest_min <- pre_process('/gpfs/commons/home/tlin/output/prs/bellenguez/updateRSID/interested_SNP/minPIP/merged_updateRSID_qc_interested_SNP.tsv')

bellenguez_qc_individual <- pre_process('/gpfs/commons/home/tlin/output/prs/bellenguez/fixed_0224/bellenguez_qc_on_individual.tsv')
bellenguez_qc_variant <- pre_process('/gpfs/commons/home/tlin/output/prs/bellenguez/fixed_0224/bellenguez_qc_on_variant.tsv')



### kunkle_qc ----
kunkle_APOE <- pre_process("/gpfs/commons/home/tlin/output/prs/kunkle/fixed_0224/APOE.tsv")
kunkle_withoutAPOE_qc <- pre_process("/gpfs/commons/home/tlin/output/prs/kunkle/fixed_0224/remove_APOE_qc_on_variant_sumstat.tsv")
kunkle_withoutAPOE <- pre_process("/gpfs/commons/home/tlin/output/prs/kunkle/fixed_0224/remove_APOE.tsv")

kunkle_cT <- pre_process('/gpfs/commons/home/tlin/output/prs/kunkle/fixed_0224/kunkle_no_qc.tsv')
kunkle_qc <- pre_process('/gpfs/commons/home/tlin/output/prs/kunkle/fixed_0224/kunkle_qc.tsv')
kunkle_qc_target <-pre_process('/gpfs/commons/home/tlin/output/prs/kunkle/fixed_0224/kunkle_qc_on_target.tsv')
kunkle_qc_base <- pre_process('/gpfs/commons/home/tlin/output/prs/kunkle/fixed_0224/kunkle_qc_on_base.tsv')

kunkle_qc_maf <- pre_process('/gpfs/commons/home/tlin/output/prs/kunkle/fixed_0224/kunkle_qc_all_maf01.tsv')
kunkle_qc_target_maf <- pre_process('/gpfs/commons/home/tlin/output/prs/kunkle/fixed_0224/kunkle_qc_target_maf01.tsv')

kunkle_qc_variant <- pre_process('/gpfs/commons/home/tlin/output/prs/kunkle/fixed_0224/kunkle_qc_on_variant_update.tsv')
kunkle_qc_individual <- pre_process('/gpfs/commons/home/tlin/output/prs/kunkle/fixed_0224/kunkle_qc_on_individual_update.tsv')
kunkle_qc_variant_sumstat <- pre_process('/gpfs/commons/home/tlin/output/prs/kunkle/fixed_0224/qc_on_variant_sumstat.tsv')
kunkle_new_beta <- pre_process('/gpfs/commons/home/tlin/output/prs/kunkle/fixed_0224/new_beta_noqc.tsv')
kunkle_new_beta <- kunkle_new_beta[-31]


### wightman ----
wightman_cT <- pre_process('/gpfs/commons/home/tlin/output/prs/wightman/before_qc.tsv')
wightman_qc <- pre_process('/gpfs/commons/home/tlin/output/prs/wightman/qc.tsv')
wightman_qc_target <- pre_process('/gpfs/commons/home/tlin/output/prs/wightman/qc_on_target.tsv')
wightman_qc_base <- pre_process('/gpfs/commons/home/tlin/output/prs/wightman/qc_on_base.tsv')
wightman_qc_individual <- pre_process('/gpfs/commons/home/tlin/output/prs/wightman/qc_on_individual_update.tsv')
wightman_qc_variant <- pre_process('/gpfs/commons/home/tlin/output/prs/wightman/qc_on_variant_update.tsv')
wightman_qc_variant_sumstat <- pre_process('/gpfs/commons/home/tlin/output/prs/wightman/qc_on_variant_sumstat.tsv')

## wightman new beta
wightman_beta <- pre_process("/gpfs/commons/home/tlin/output/prs/wightman/new_beta_max_snp_10_noclump_pT.tsv")
wightman_beta_cT <- pre_process("/gpfs/commons/home/tlin/output/prs/wightman/new_beta_max_snp_10_clump_pT.tsv")
#wightman_beta_polypred <- pre_process("/gpfs/commons/home/tlin/output/prs/wightman/susie_new_beta_max_snp_10_polypred.tsv")  ## only the last step is using polypred
wightman_susie_max10_polypred<- pre_process('/gpfs/commons/home/tlin/output/prs/wightman/susie_new_beta_max_snp_10_polypred.tsv')


###  polyfun-Pred ----
kunkle_polypred <- pre_process('/gpfs/commons/home/tlin/output/prs/polypred/kunkle/kunkle_polypred.tsv')
bellenguez_polypred_new <- pre_process('/gpfs/commons/home/tlin/output/prs/polypred/bellenguez/fixed_0224_polypred.prs')
wightman_polypred <- pre_process('/gpfs/commons/home/tlin/output/prs/polypred/wightman/fixed_0224.prs.tsv')

## new plink ------
kunkle_adsp <- pre_process('/gpfs/commons/home/tlin/output/prs/new_plink/kunkle/kunkle_ADSP.tsv')
bellenguez_adsp <- pre_process('/gpfs/commons/home/tlin/output/prs/new_plink/bellenguez/bellenguez_ADSP.tsv')
wightman_adsp <- pre_process('/gpfs/commons/home/tlin/output/prs/new_plink/wightman/wightman_ADSP.tsv')


kunkle_adsp_qc <- pre_process('/gpfs/commons/home/tlin/output/prs/new_plink/kunkle/kunkle_ADSP_qc_all.tsv')
bellenguez_adsp_qc <- pre_process('/gpfs/commons/home/tlin/output/prs/new_plink/bellenguez/bellenguez_ADSP_qc_all.tsv')
wightman_adsp_qc <- pre_process('/gpfs/commons/home/tlin/output/prs/new_plink/wightman/wightman_ADSP_qc_all.tsv')

kunkle_UKBB_qc <- pre_process('/gpfs/commons/home/tlin/output/prs/new_plink/kunkle/kunkle_ADSP_UKBB_qc.tsv')
bellenguez_UKBB_qc <- pre_process('/gpfs/commons/home/tlin/output/prs/new_plink/bellenguez/bellenguez_ADSP_UKBB_qc.tsv')
wightman_UKBB_qc <- pre_process('/gpfs/commons/home/tlin/output/prs/new_plink/wightman/wightman_ADSP_UKBB_qc.tsv')



### susie -----
kunkle_susie <- pre_process("/gpfs/commons/home/tlin/output/prs/polypred/kunkle/susie.prs.tsv")
susie <- rename_preprocess("/gpfs/commons/home/tlin/output/prs/bellenguez_susie/bellenguez_susie_prs_PC.tsv")

## Others
PRSice <- rename_preprocess("/gpfs/commons/home/tlin/output/prs/PRSice_pheno.tsv")
sbayesR = rename_preprocess("/gpfs/commons/home/tlin/output/prs/sbayesR.tsv")

##density_plot -------
density_plot(sbayesR, name = "sbayesR")
par(mfrow=c(1,3))

p_threshold_col <- c("PRS_005","PRS_05","PRS_1","PRS_5")
pT_subset = pT[p_threshold_col]

for(i in c(33,35,37,34,36,38)) {       # for-loop over columns
  plot(density((pT[, i])),col='grey', lty = 2, lwd = 2,main=str_replace(colnames(pT)[i],"PRS_", "p=0."), xlab='PRS')
  lines(density(pT[pT$Diagnosis == 1,i]), col = 'red')
  lines(density(pT[pT$Diagnosis == 0,i]), col = "blue")
  legend("topright", legend=c("ALL","Case", "Control"),
         col=c("darkgrey","red", "blue"), lty=c(2,1,1), cex=0.8,
         box.lty=0)
}

mtext("bellenguez (updated RSID)",                   # Add main title
      side = 3,
      line = -1.25,
      outer = TRUE)


for(i in c(3,5,7)) {       # for-loop over columns
  plot(density((kunkle_cT[, i])),col='grey', lty = 2, lwd = 2,main=str_replace(colnames(kunkle_cT)[i],"PRS_", "p=0."), xlab='PRS')
  lines(density(kunkle_cT[kunkle_cT$Diagnosis == 1,i]), col = 'red')
  lines(density(kunkle_cT[kunkle_cT$Diagnosis == 0,i]), col = "blue")
  legend("topright", legend=c("ALL","Case", "Control"),
         col=c("darkgrey","red", "blue"), lty=c(2,1,1), cex=0.8,
         box.lty=0)
}
mtext("kunkle",                   # Add main title
      side = 3,
      line = -1.25,
      outer = TRUE)



##draw max_snp_per ld
prs_col =  c("PRS1","PRS3","PRS5","PRS10")

##plotting case/control plot ---- 
## However, instead of plotting it in whole population, we should plot them with race_separate plots (facet plot)
## Hence, this function is not useful compare with the facet one.
plot_density_diag <- function(df, col, title){
  for(i in 1:length(col)) {  
    if(grepl("_", col[1])){
      if(grepl("e5", col[i]))
        plot(density((df[, col[i]])),col="grey", lty = 2, lwd = 2, main=str_replace(col[i],"PRS_", "pT = "), xlab='PRS')
      else
        plot(density((df[, col[i]])),col="grey", lty = 2, lwd = 2, main=str_replace(col[i],"PRS_", "pT = 0."), xlab='PRS')
    }
    else
      plot(density((df[, col[i]])),col="grey", lty = 2, lwd = 2, main=str_replace(col[i],"PRS", "max_snp_per_LD = "), xlab='PRS')
    lines(density(df[df$Diagnosis == 1,col[i]]), col = 'red')
    lines(density(df[df$Diagnosis == 0, col[i]]), col = 'blue')
    legend("topright", legend=c("ALL","Case", "Control"),
           col=c("darkgrey","red", "blue"), lty=c(2,1,1), cex=0.8,
           box.lty=0)
  }
  mtext(title,                
        side = 3,line = -1.25,outer = TRUE)
}
par(mfrow=c(4,1))
plot_density_diag(kunkle_cT, c("PRS_e5", "PRS_001", "PRS_005", "PRS_01", "PRS_05", "PRS_1", "PRS_5"),"kunkle_cT")
plot_density_diag(susie,c("PRS1","PRS10"), "susie")
plot_density_diag(kunkle_susie,c("PRS1","PRS3","PRS7","PRS10"), "kunkle_susie")
plot_density_diag(kunkle_new_beta,c("PRS_e5", "PRS_005", "PRS_05", "PRS_5"), "kunkle_pT_newbeta")

## create a race-separated pivot df for poly_pred (which is easier for plotting). Where method = how many causal SNPs are assumed in the LD block.
pivot_df<- function(df){ 
  df_new <- df %>%
    subset(final_population ==  "EUR" | final_population == "AMR"| final_population =="AFR") %>%   ## only keep the race of interest
    pivot_longer(                           ## can be apply to different methods (polyfunpred or pT)
      cols = starts_with("PRS"),
      names_to = 'method', 
      values_to = 'PRS') %>% 
    as.data.frame()
  #df_new$Diagnosis = as.factor(df_new$Diagnosis)
  df_new$Diagnosis_state=ifelse(df_new$Diagnosis == 0, "control","case")  ## create a new column so that it will show case/control when plotting
  return(df_new)
}
plot_facut_dist <- function(df, title){
  ggplot(pivot_df(df), aes(x =PRS, group=Diagnosis_state,fill = Diagnosis_state, color = Diagnosis_state))+
    geom_density(alpha = 0.01) +
    facet_grid(method~final_population)+theme_bw()+ggtitle(title)
} ## plot race/method separated dist plots

plot_facut_dist(kunkle_polypred, "kunkle_polypred")
plot_facut_dist(kunkle_new_beta, "kunkle_new_beta")
plot_facut_dist(kunkle_susie, "kunkle_susie")

plot_facut_dist(wightman_beta, "wightman_new_beta")
plot_facut_dist(wightman_beta_cT, "wightman_new_beta_cT")
plot_facut_dist(wightman_susie_max10_polypred, "wightman_susie_new_beta")


## plot it another way? -----
diagnosis_label = c("Control","Case")
names(diagnosis_label) = c(0,1)

ggplot(test, aes(x = PRS, group = pT,fill = pT, color =pT))+
  geom_density(alpha = 0.3) +
  facet_grid(Diagnosis~final_population, 
             labeller = labeller(Diagnosis= diagnosis_label))+
  theme_bw()+ggtitle("kunkle_cT")

### plot all of each method----
par(mfrow=c(1,2))

plot(density(pT$PRS_5),col="burlywood3", lty = 2, lwd = 2, xlab='PRS', main = " ", xlim = c(-0.016,0.004))
lines(density(bellenguez_update$PRS10), col="darkolivegreen3",lwd = 2, lty = 2)
abline(v = mean(pT$PRS_5), col = "burlywood3")
abline(v = mean(bellenguez_update$PRS10), col = "darkolivegreen3")

lines(density(PRSice$PRS), col = 'cadetblue1', lwd =2, lty = 2)
abline(v = mean(PRSice$PRS), col = "cadetblue1")

text(mean(pT$PRS_5)+0.002,1100,round(mean(pT$PRS_5),4), col = 'brown')
text(mean(pT$PRS_5)+0.004,1200,'clumping+ pT (p = 0.5)', col = 'brown')
text(mean(bellenguez_update$PRS10),400,round(mean(bellenguez_update$PRS10),4), col = 'darkolivegreen4')
text(mean(bellenguez_update$PRS10),500,"polypred (max_snp = 10)", col = 'darkolivegreen4')
text(mean(PRSice$PRS)-0.002,1400,"PRSice", col = 'deepskyblue')
text(mean(PRSice$PRS)-0.002,1300,round(mean(PRSice$PRS),6), col = 'deepskyblue')

## SbayesR & PRSice
plot(density(sbayesR$PRS), col = 'darkseagreen3', lwd =2 ,lty = 2,xlab = 'PRS', xlim=c(-1e-03,6e-05), main = ' ')
abline(v = mean(sbayesR$PRS), col = "darkseagreen3")
lines(density(PRSice$PRS), col = 'cadetblue1', lwd =2, lty = 2)
abline(v = mean(PRSice$PRS), col = "cadetblue1")

text(mean(PRSice$PRS),2200,"PRSice", col = 'deepskyblue')
text(mean(PRSice$PRS),1000,round(mean(PRSice$PRS),6), col = 'deepskyblue')

text(mean(sbayesR$PRS),20000,"SBayesR", col = 'darkcyan')
text(mean(sbayesR$PRS),19000,round(mean(sbayesR$PRS),6), col = 'darkcyan')

  
## Extract specific race 
extract_eur <-function (df){
  EUR = df[df$final_population == "EUR",]
  return(EUR)
}
extract_afr <- function (df){
  AFR = df[df$final_population == "AFR",]
  return(AFR)
}
extract_amr <- function (df){
  AMR = df[df$final_population == "AMR",]
  return(AMR)
}


col_roc <- list("PRS_001","PRS_005","PRS_01","PRS_05","PRS_1","PRS_5")
col_roc_E5 <- list("PRS_e5","PRS_001","PRS_005","PRS_01","PRS_05","PRS_1","PRS_5")

##AUC function ------
plot_roc <- function(roclist, roccol,legend="PRS_", replace="p = 0.",title = FALSE){
  legendname = list()
  for (i in roccol){
    if(i == "PRS_e5"){
      legendname = append(legendname,paste('p = 1e-5, auc =',round(roclist[[i]]$auc,3)))
    }else{
      legendname = append(legendname,paste(str_replace(i,legend,replace),', auc =',round(roclist[[i]]$auc,3)))
    }
  }
  print(legendname)
  curve <- ggroc(roclist)
  curve + ggtitle(title)+theme_bw()+
    scale_colour_discrete(name="PRS",labels = c(legendname))
}   
roc_result <-function(df, title=' ',column_for_roc, plot=TRUE){
  if('PRS_e5' %in% column_for_roc){
    roc_list <- roc(Diagnosis ~ PRS_e5+PRS_001+PRS_005+PRS_01+PRS_05+PRS_1+PRS_5, data = df)
  }else if ('PRS10' %in% column_for_roc){
    roc_list <- roc(Diagnosis ~ PRS1+PRS3+PRS5+PRS7+PRS10, data = df)
  }else if ('new_beta' %in% column_for_roc){
    roc_list <- roc(Diagnosis ~ new_beta + new_beta_cT +new_beta_susie, data=df)
  }else{
    roc_list <- roc(Diagnosis ~ PRS_001+PRS_005+PRS_01+PRS_05+PRS_1+PRS_5, data = df)
  }
  if(plot == TRUE){
    plot_roc(roc_list,column_for_roc, title=title)
  }
  else{
    len = length(roc_list)
    my_list=list()
    for (i in 1:len){
      my_list[i] = roc_list[[i]]$auc
    }
    return(my_list)
  }
}

plot_ethnic_roc <- function(df, title, col, plot=TRUE){
  #all = roc_result(df, title = paste(title, ", all population") , column_for_roc = col, plot=plot)
  EUR = roc_result(extract_eur(df), title = paste(title, ", EUR") , column_for_roc = col,plot=plot)
  AFR = roc_result(extract_afr(df), title = paste(title, ", AFR") , column_for_roc = col, plot=plot)
  AMR = roc_result(extract_amr(df), title = paste(title, ", AMR") , column_for_roc = col, plot=plot)
  
  if(plot == TRUE){
    plot_grid(EUR, AFR, AMR,ncol = 1, nrow = 3)
  }
  else{
    df <- data_frame(matrix(ncol = 0, nrow = 7))
    df$PRS = unlist(col)
    df$EUR = unlist(EUR)
    df$AFR = unlist(AFR)
    df$AMR = unlist(AMR)
    return(df[,c(2:5)])
  }
}


plot_auc_cT <- function(df, df_qc, df_qc_base, df_qc_target, sumstat_name, col_to_use, eur=FALSE){
  
  if(eur==TRUE){
    df = extract_eur(df)
    df_qc = extract_eur(df_qc)
    df_qc_base = extract_eur(df_qc_base)
    df_qc_target = extract_eur(df_qc_target)
  }
    
    
  cT_plot = roc_result(df, title = paste(sumstat_name, ", c+pT") , column_for_roc = col_to_use)
  qc_plot = roc_result(df_qc, title = paste(sumstat_name, ", QC") , column_for_roc = col_to_use)
  qc_on_base_plot = roc_result(df_qc_base, title = paste(sumstat_name, ", QC on summary stat") , column_for_roc = col_to_use)
  qc_on_target_plot = roc_result(df_qc_target, title = paste(sumstat_name, "QC on target") , column_for_roc = col_to_use)
  
  plot_grid(cT_plot, qc_plot, qc_on_base_plot, qc_on_target_plot,ncol = 2, nrow = 2)
  
}

# c+T (before and after qc)-----
## bellenguez-----

##bellenguez_qc
plot_auc_cT(bellenguez_cT, bellenguez_qc, bellenguez_qc_on_base, bellenguez_qc_on_target, "bellenguez(updateRSID)", col_roc)
plot_auc_cT(bellenguez_cT, bellenguez_qc, bellenguez_qc_on_base, bellenguez_qc_on_target, "bellenguez(updateRSID, EUR)", col_roc, eur =TRUE)


## bellenguez_ethnics ----
plot_ethnic_roc(bellenguez_cT, title= 'bellenguez(updateRSID), pT', col_roc)
plot_ethnic_roc(bellenguez_adsp, title= 'bellenguez(new_plink), pT', col_roc)
plot_ethnic_roc(bellenguez_adsp_qc, title= 'bellenguez(new_plink, qc), pT', col_roc)


##bellenguez_interested_SNP ------
auc(roc(Diagnosis~PRS, data = bellenguez_interested))  # 0.5187
auc(roc(Diagnosis~PRS, data = bellenguez_interest_max))  # 0.515
auc(roc(Diagnosis~PRS, data = bellenguez_interest_min))  #0.5109

auc(roc(Diagnosis~PRS, data = extract_eur(bellenguez_interested))) ##0.5308
auc(roc(Diagnosis~PRS, data = extract_eur(bellenguez_interest_max))) ##0.5247
auc(roc(Diagnosis~PRS, data = extract_eur(bellenguez_interest_min))) ##0.5224

auc(roc(Diagnosis~PRS, data = extract_afr(bellenguez_interested))) #0.529
auc(roc(Diagnosis~PRS, data = extract_amr(bellenguez_interested))) #0.4991

auc(roc(Diagnosis~PRS, data = bellenguez_qc_interested)) # 0.5255
auc(roc(Diagnosis~PRS, data = extract_eur(bellenguez_qc_interested))) ##0.5422
auc(roc(Diagnosis~PRS, data = extract_afr(bellenguez_qc_interested))) #0.5215
auc(roc(Diagnosis~PRS, data = extract_amr(bellenguez_qc_interested))) #0.5182

auc(roc(Diagnosis~PRS10, data = susie)) ##0.5075
auc(roc(Diagnosis~PRS10, data = bellenguez_fixed_0224)) #0.5223

auc(roc(Diagnosis~PRS1, data = kunkle_susie)) ## 0.5046
auc(roc(Diagnosis~PRS1, data = kunkle_polypred)) # 0.5026


## kunkle------

## Kunkle
plot_auc_cT(kunkle_cT, kunkle_qc, kunkle_qc_base, kunkle_qc_target, "kunkle", col_roc_E5)
plot_auc_cT(kunkle_cT, kunkle_qc, kunkle_qc_base, kunkle_qc_target, "kunkle, EUR", col_roc_E5, eur=TRUE)

### Qced

NoAPOE_QC = roc_result(kunkle_withoutAPOE_qc,title="kunkle without APOE, qc on sumstats and variants", column_for_roc=col_roc_E5)
NoAPOE = roc_result(kunkle_withoutAPOE,title="kunkle without APOE", column_for_roc=col_roc_E5)

plot_grid(kunkle_pt_plot,kunkle_pt_qc_plot_variant_sumstat,NoAPOE,NoAPOE_QC,ncol = 2, nrow = 2)

kunkle_pt_plot=roc_result(kunkle_cT,title="kunkle (no qc)", column_for_roc = col_roc_E5)
kunkle_pt_qc_plot_target=roc_result(kunkle_qc_target_maf,title="QC on target", column_for_roc = col_roc_E5)
kunkle_pt_qc_plot_variant=roc_result(kunkle_qc_variant,title="QC on variant", column_for_roc = col_roc_E5)
kunkle_pt_qc_plot_variant_sumstat=roc_result(kunkle_qc_variant_sumstat,title="QC on variant and sumstat", column_for_roc = col_roc_E5)
plot_grid(kunkle_pt_qc_plot,kunkle_pt_qc_plot_target, kunkle_pt_qc_plot_variant,kunkle_pt_qc_plot_variant_sumstat, ncol = 2, nrow = 2)

kunkle_pt_qc_plot=roc_result(extract_eur(kunkle_cT),title="kunkle, EUR, (no qc)", column_for_roc = col_roc_E5)
kunkle_pt_qc_plot_target=roc_result(extract_eur(kunkle_qc_target_maf),title="EUR, QC on target", column_for_roc = col_roc_E5)
kunkle_pt_qc_plot_variant=roc_result(extract_eur(kunkle_qc_variant),title="EUR, QC on variant", column_for_roc = col_roc_E5)
kunkle_pt_qc_plot_variant_sumstat=roc_result(extract_eur(kunkle_qc_variant_sumstat),title="EUR, QC on variant and sumstat", column_for_roc = col_roc_E5)
plot_grid(kunkle_pt_qc_plot,kunkle_pt_qc_plot_target, kunkle_pt_qc_plot_variant,kunkle_pt_qc_plot_variant_sumstat, ncol = 2, nrow = 2)

kunkle_pt_qc_plot_maf=roc_result(kunkle_qc_maf,title="QC, MAF >0 0.1%", column_for_roc = col_roc_E5)
##kunkle_pt_qc_plot_base=roc_result(kunkle_qc_base, title="kunkle_qc_summary_stats",column_for_roc = col_roc_E5)
plot_grid(kunkle_pt_qc_plot,kunkle_pt_qc_plot_maf,ncol = 2, nrow = 2)



##cros group -----
par(mfrow=c(3,1))
col_roc_polypred <- list("PRS1","PRS3","PRS5","PRS7","PRS10")
plot_ethnic_roc (kunkle_cT, 'kunkle', col_roc_E5)
plot_ethnic_roc(kunkle_qc_variant_sumstat, "kunkle_qc_sumstat_variant", col_roc_E5)
plot_ethnic_roc (kunkle_polypred, 'kunkle', col_roc_polypred)
plot_ethnic_roc (wightman_polypred, 'wightman', col_roc_polypred)
plot_ethnic_roc (bellenguez_fixed_0224,'bellenguez', col_roc_polypred)
auc(roc(Diagnosis~PRS, data = extract_eur(bellenguez_polypred_new)))
auc(roc(Diagnosis~PRS, data = extract_afr(bellenguez_polypred_new)))
auc(roc(Diagnosis~PRS, data = extract_amr(bellenguez_polypred_new)))
plot_ethnic_roc(kunkle_new_beta, 'kunkle, (effectsize * pip) ',col_roc_E5)
plot_ethnic_roc(kunkle_susie, "kunkle(susie)", col_roc_polypred)
plot_ethnic_roc(susie, "bellenguez(susie)", col_roc_polypred)

new_beta_col = list("new_beta","new_beta_cT","new_beta_susie")
plot_ethnic_roc(wightman_new_beta, 'wightman, newbeta(effectsize * pip) ',new_beta_col, plot=TRUE)
title='wightman, newbeta(effectsize * pip) '

roc(Diagnosis ~ PRS_new_beta+new_beta_cT+new_beta_susie, df=wightman_new_beta)

roc(Diagnosis~PRS_new_beta+new_beta_cT+new_beta_susie, data = extract_eur(wightman_new_beta))


wightman_new_beta=wightman_beta
wightman_new_beta["new_beta_cT"] = wightman_beta_cT$PRS
wightman_new_beta["new_beta_susie"] = wightman_susie_max10_polypred$PRS
wightman_new_beta$PRS
colnames(wightman_new_beta)[which(names(wightman_new_beta) == 'PRS')] <- 'new_beta'
colnames(wightman_new_beta)[which(names(wightman_new_beta) == 'PRS_new_beta')] <- 'new_beta'
PRS_new_beta
wightman_new_beta

## APOE
plot_ethnic_roc(kunkle_withoutAPOE, "kunkle remove APOE region", col_roc_E5)
col_roc <- list("only2SNP","no2SNP","no2SNP_qc")
kunkle_APOE_roc <- roc(Diagnosis ~ only2SNP+no2SNP+no2SNP_qc, data = kunkle_APOE)
auc(roc(Diagnosis~only2SNP, data = kunkle_APOE)) ##0.621
auc(roc(Diagnosis~no2SNP, data = kunkle_APOE))##0.5131

auc(roc(Diagnosis~only2SNP, data = extract_eur(kunkle_APOE))) ##0.6187
auc(roc(Diagnosis~no2SNP, data = extract_eur(kunkle_APOE)))  ##0.5834
auc(roc(Diagnosis~no2SNP, data = extract_amr(kunkle_APOE)))
auc(roc(Diagnosis~no2SNP, data = extract_afr(kunkle_APOE)))

plot_ethnic_roc(kunkle_withoutAPOE_qc, 'kunkle (No APOE, qc on sumstats and variants) ',col_roc_E5)
plot_ethnic_roc(kunkle_withoutAPOE, 'kunkle (No APOE)',col_roc_E5)

## wightman ----
plot_auc_cT(wightman_cT,wightman_qc, wightman_qc_base, wightman_qc_target, 'Wightman', col_roc_E5)
plot_ethnic_roc(wightman_cT, "wightman_qc", col_roc_E5)


## polyfun----

##susie

susie_roc = plot_roc(roc_list,col_roc, legend = 'PRS', replace = 'Max SNP = ', title="bellenguez(SuSiE)")
plot_grid(bellenguez_roc,susie_roc,ncol = 2, nrow = 1)
plot_grid(bellenguez_roc,pt_plot,ncol= 2, nrow = 1)

## all
par(mfrow=c(1,1))
color = c("red","blue","darkgreen","slategray")
plot(roc(sbayesR$Diagnosis~sbayesR$PRS), col = color[1],main="auROC")
lines(roc(PRSice$Diagnosis~PRSice$PRS),col =  color[2])
lines(roc(bellenguez_update$Diagnosis~bellenguez_update$PRS5),col = color[3])
lines(roc(pT$Diagnosis~pT$PRS_5),col = color[4])


legend("bottom", legend=c(paste("SbayesR, auc = ",round(roc(sbayesR$Diagnosis~sbayesR$PRS)$auc,4)),
                               paste("PRSice, auc = ",round(roc(PRSice$Diagnosis~PRSice$PRS)$auc,4)),
                               paste("Polyfun-Pred, auc = ",round(roc(bellenguez_update$Diagnosis~bellenguez_update$PRS5)$auc,4)),
                               paste("Clumping+pT, auc = ",round(roc(pT$Diagnosis~pT$PRS_5)$auc,4))),
                      col=color, lty=1, cex=0.8,box.lty=0)
     

## Regression function----
## plot beta and pvalue
p_beta_plot <- function(mod, prs_pos, title){
  pvalue = ggplot(mod[prs_pos,], aes(x =rownames(mod[prs_pos,]), y =  Pr...z..)) + ylab('P value')+xlab(' ')+
    ggtitle(title)+
    geom_point( size=2, shape=17,fill="white")+ylim(c(min(mod[,3])-0.1,max(mod[,3])+0.1))+
    scale_fill_brewer(palette="Paired")+
    geom_text(aes(label=Pr...z..), vjust=1.6, color="black",
              position = position_dodge(1), size=3)+  theme_minimal()
  
  beta = ggplot(mod[prs_pos,], aes(x =rownames(mod[prs_pos,]), y =  Estimate, fill = PRS)) + ylab('Beta')+xlab(' ')+
    geom_bar(stat="identity") +scale_fill_brewer(palette="Paired")+
    geom_errorbar(aes(ymin=Estimate-Std..Error, ymax=Estimate+Std..Error), width=.2,
                  position=position_dodge(.9),) +
    geom_text(aes(label=Estimate), vjust=1.2, color="darkred",
              position = position_dodge(1), size=3)+  theme_minimal()
  plot_grid(pvalue, beta)
}

##calcualte pseudo rsquare
rsq <- function(formula, data, indices=FALSE) {    
  if (indices !=FALSE){
    d <- data[indices,] # allows boot to select sample
  }
  else 
    d = data 
  fit <- glm(formula = formula,family = 'binomial', data = d)
  return(RsqGLM(fit, plot = FALSE)$Nagelkerke)
}

##resample
log_reg <- function(df,prs,main_title, plot = TRUE, legend="PRS_", boot_num = FALSE, replace="pT = 0."){
  set.seed(10)
  mod1 <- glm(Diagnosis ~ Sex + Age, data= df, family=binomial)
  mod1_R2 <- RsqGLM(mod1, plot=FALSE)$Nagelkerke
  for (i in 1:length(prs)){
    frm <- as.formula(paste("Diagnosis ~ Sex + Age+", prs[[i]]))
    if(boot_num != FALSE){
      mod <- boot(data=df, statistic=rsq,
                     R=boot_num, formula=frm)
      #CI  <- boot.ci(boot.out=mod,type='norm')$normal[c(2,3)]
     
      sd = sd(mod$t)
      mean =  mean(mod$t)
      CI = c(mean-2*sd, mean+2*sd)
      ##table
      if(i ==1){
        R2 <- data.frame (
                          'PRS' = prs[[i]],
                          'boot_mean'  = mean,
                          'boot_CI_upper' = CI[2],
                          'boot_CI_lower' = CI[1],
                          'SD' = sd,
                          "null_R2" = mod1_R2)
      }
      else{
        R2_else <- data.frame ('PRS' = prs[[i]],
                          'boot_mean'  = mean,
                          'boot_CI_upper' = CI[2],
                          'boot_CI_lower' = CI[1],
                          'SD' = sd,
                          "null_R2" = mod1_R2)
        R2 <- rbind(R2, R2_else)
        if(dim(R2)[1]==length(prs)){
          if(plot == TRUE)
            plotR2_boot(R2,main_title, prs,boot_num)
          return (R2)
        }
      }
    }  ##bootstrap or not
    else{
      mod <- glm(formula = frm,family = 'binomial', data = df)
      mod_R2 <-  RsqGLM(mod, plot=FALSE)$Nagelkerke 
      ## table
      if(i ==1){  ## if theres only one threshold
        cof <- data.frame(round(summary(mod)$coefficients[,c(1,2,4)],4))
        cof['null_R2'] = round(mod1_R2,5)
        cof['R2_prs'] = round(mod_R2,5)
        cof['PRS'] = prs[i]
      }
      else{
        cof2 <- data.frame(round(summary(mod)$coefficients[,c(1,2,4)],4))
        cof2['null_R2'] = round(mod1_R2,5)
        cof2['R2_prs'] = round(mod_R2,5)
        cof2['PRS'] =prs[i]
        cof <- rbind(cof, cof2)
        #print(cof)
        if(dim(cof)[1]/4==length(prs) ){
          if (plot == TRUE)
            plotR2_boot(cof,main_title, prs,boot_num)
          return (cof)
        }
      }
    }
  }
}


log_reg(kunkle_cT,col_roc, " ", FALSE, boot_num=5)
par(xpd = FALSE,mfrow=c(1,1))
plotR2_boot <- function(log_output, header, prs_col, boot){
  par(xpd = FALSE)
  if(boot == FALSE){   ## no bootstrapping
    prs_col = seq(from=1, to=dim(log_output)[1], by=4)
    plot(log_output$R2_prs[prs_col], ylim  = c(min(log_output$null_R2[1], min(log_output$R2_prs[prs_col]))*0.95, max(log_output$null_R2[1], max(log_output$R2_prs[prs_col]))),
         ylab='R-squared', xlab = "different PRS",pch=4,col="darkgreen", main=header,xaxt = "n")
    text(seq(from=1,to=length(prs_col)),log_output$R2_prs[prs_col] , 
         labels = log_output$R2_prs[prs_col], adj = c(0.5,2), xpd = TRUE, cex=1, col="darkgreen") 
    axis(1,                         # Define x-axis manually
         at = 1:length(prs_col),
         labels = c(str_replace(log_output$PRS[1],"PRS_e5","p = 1*e-5"),str_replace(log_output$PRS[seq(from=1, to=dim(log_output)[1], by=4)][-1],"PRS_","p = 0.")))
    
    
  }
  else
  {
    plot(log_output$boot_mean, ylim  = c(min(log_output$null_R2[1], min(log_output$boot_mean))*0.1, max(log_output$boot_mean[1], max(log_output$boot_mean)*1.5)),
         ylab='R-squared', xlab = "different PRS",pch=4,col="darkgreen", main=header,xaxt = "n")
    arrows(1:7, log_output$boot_CI_upper, 1:7, log_output$boot_CI_lower, length=0.05, angle=90, code=3, col='grey58') 
    text(seq(from=1,to=length(prs_col)),log_output$boot_mean , 
         labels = round(log_output$boot_mean,5), adj = c(0.5,2), xpd = TRUE, cex=1, col="darkgreen") 
    axis(1,                         # Define x-axis manually
         at = 1:length(prs_col),
         labels = c(str_replace(prs_col[1],"PRS_e5","p = 1*e-5"),str_replace(prs_col[seq(from=1, to=length(prs_col))][-1],"PRS_","p = 0.")))
  }
  
  abline(h=log_output$null_R2[1], col=c("blue"), lty=2, lwd=2)
  text((length(prs_col)+1)/2,log_output$null_R2[1] , labels =  paste("Null model: ",round(log_output$null_R2[1],5)), adj = c(0.5,-1), xpd = TRUE, cex=1, col="blue")
  
}

plot_ethnic_R2 <- function(df, col, title, boot_num, replace='', plot= TRUE){
  #all = log_reg(df, col, paste(title, ", ALL"), plot =T, legend="PRS",boot_num = boot_num, replace=replace)
  print(paste("running EUR with bootstrapping ", boot_num, ' times'))
  EUR = log_reg(extract_eur(df), col,  paste(title, ", EUR"), plot = plot, legend="PRS",boot_num = boot_num, replace=replace)
  EUR$ethnicity="EUR"
  print(paste("running AFR with bootstrapping ", boot_num, ' times'))
  AFR = log_reg(extract_afr(df), col, paste(title, ", AFR"), plot = plot, legend="PRS",boot_num = boot_num, replace=replace)
  AFR$ethnicity="AFR"
  print(paste("running AMR with bootstrapping ", boot_num, ' times'))
  AMR = log_reg(extract_amr(df), col,  paste(title, ", AMR"), plot = plot, legend="PRS",boot_num = boot_num, replace=replace)
  AMR$ethnicity="AMR"
  if(plot != "TRUE"){
    return(rbind(EUR, AFR, AMR))
  }
  
}


# test pseudo rsquare
NagelkerkeR2(mod2)$R2
RsqGLM(mod2)$Nagelkerke
# test ti get 95% confidence interval
boot.ci(results, type="bca")
result = boot.ci(boot.out=results,type='norm')
print(c(mean(results$t) - 2*sd(results$t), mean(results$t) + 2*sd(results$t)))

## wrapped up function ## for easier calculation for pseudo r2
plot_R2_cT <- function(df, df_qc, df_qc_base, df_qc_target, sumstat_name, col_to_use, plot=TRUE, boot_num = FALSE, eur =FALSE){
  par(mfrow=c(2,2),xpd=FALSE)
  if (eur ==TRUE){
    df = extract_eur(df)
    df_qc = extract_eur(df_qc)
    df_qc_base = extract_eur(df_qc_base)
    df_qc_target = extract_eur(df_qc_target)
  }

  cT= log_reg(df, col_to_use, paste(sumstat_name, ", c+pT") , plot = plot,  boot_num = boot_num,replace="max_snp_per_locus=")
  qc= log_reg(df_qc, col_to_use, paste(sumstat_name, ", QC") , plot = plot,  boot_num = boot_num,replace="max_snp_per_locus=")
  qc_on_base= log_reg(df_qc_base, col_to_use, paste(sumstat_name, ", QC on summary stat") , plot = plot,  boot_num = boot_num,replace="max_snp_per_locus=")
  qc_on_target= log_reg(df_qc_target, col_to_use, paste(sumstat_name, ", QC on target") , plot = plot,  boot_num = boot_num,replace="max_snp_per_locus=")
}


##plot facet
##R2
## plot ethnic facut plot for one sumstat
plot_ethnic_R2_facut <- function(QC1, QC2, QC3, col, boot_num, title,QC1name="QC_on_base", QC2name="QC_on_base+variants",QC3name="no_qc"){
  QC1 = plot_ethnic_R2(QC1, col, '', boot_num, plot="F")
  QC1$qc_status = QC1name
  QC2 = plot_ethnic_R2(QC2, col, '', boot_num, plot="F")
  QC2$qc_status = QC2name
  QC3 = plot_ethnic_R2(QC3, col, '', boot_num, plot="F")
  QC3$qc_status = QC3name
  all = rbind(QC1,QC2, QC3)
  all %>% 
    filter(PRS !="PRS-001" & PRS != "PRS_1") 
  
  all$PRS = str_replace(all$PRS, 'PRS_', 'p = 0.')
  all[all$PRS=='p = 0.e5',]$PRS =" p = 1e-5"
 
  plot <- ggplot(data = all, aes(x=boot_mean, y = PRS, shape=qc_status, color = qc_status))+
    geom_point(size=3)+
    geom_pointrange(aes(xmin=boot_CI_lower, xmax=boot_CI_upper), linetype="dotted") +
    facet_wrap(~ethnicity, ncol=1)+
    xlab('R squared')+ ggtitle(title)+
    theme_bw()
  return(plot)
}

## plot ethnic facut plot for multiple sumstat
plot_R2_facut_allsumstat<- function(s1_qc1, s1_qc2, s1_qc3,s2_qc1, s2_qc2, s2_qc3,s3_qc1, s3_qc2, s3_qc3,col = col_roc_E5, QC1name="QC_on_base", QC2name="QC_on_base+variants",QC3name="no_qc",boot_num=50){
  a = plot_ethnic_R2_facut(s1_qc1, s1_qc2, s1_qc3, col, boot_num, 'Kunkle', QC1name, QC2name,QC3name)
  b = plot_ethnic_R2_facut(s2_qc1, s2_qc2, s2_qc3, col, boot_num, 'Bellenguez',QC1name, QC2name,QC3name)
  c = plot_ethnic_R2_facut(s3_qc1, s3_qc2, s3_qc3, col, boot_num,'Wightman', QC1name, QC2name,QC3name)
  prow <- plot_grid(a+ theme(legend.position="none"), 
                    b+ theme(legend.position="none",axis.text.y = element_blank())+ ylab(NULL), 
                    c+ theme(legend.position="none",axis.text.y = element_blank())+ylab(NULL),
                    ncol = 3, nrow = 1)
  legend_b <- get_legend(
    c+ guides(color = guide_legend(nrow = 1)) +
      theme(legend.position = "bottom")
  )
  plot_grid(prow, legend_b, ncol=1,rel_heights=c(3,.4))

}

## auc ----debug ------
plot_ethnic_roc_facut <- function(QC1, QC2, QC3, col, title,QC1name="QC_on_base", QC2name="QC_on_base+variants",QC3name="no_qc"){
  QC1 = plot_ethnic_roc(QC1, '',col, plot="F")
  QC1$qc_status = QC1name
  QC2 = plot_ethnic_roc(QC2, '', col, plot="F")
  QC2$qc_status = QC2name
  QC3 = plot_ethnic_roc(QC3, '', col, plot="F")
  QC3$qc_status = QC3name
  all = rbind(QC1,QC2, QC3)
  all %>% 
    filter(PRS !="PRS-001" & PRS != "PRS_1") 
  
  all$PRS = str_replace(all$PRS, 'PRS_', 'p = 0.')
  all[all$PRS=='p = 0.e5',]$PRS =" p = 1e-5"
  print(colnames(all))
  
  plot <- ggplot(data = all, aes(x=auc, y = PRS, shape=qc_status, color = qc_status))+
    geom_point(size=3)+
    facet_wrap(~ethnicity, ncol=1)+
    xlab('AUC')+ ggtitle(title)+
    theme_bw()
  return(plot)
}
plot_auc_facut_all_sumstat <- function(s1_qc1, s1_qc2, s1_qc3,s2_qc1, s2_qc2, s2_qc3,s3_qc1, s3_qc2, s3_qc3, col = col_roc_E5, QC1name="QC_on_base", QC2name="QC_on_base+variants",QC3name="no_qc"){
  a = plot_ethnic_roc_facut(s1_qc1, s1_qc2, s1_qc3, col,'Kunkle',QC1name, QC2name,QC3name)
  b = plot_ethnic_roc_facut(s2_qc1, s2_qc2, s2_qc3, col,'Bellenguez',QC1name, QC2name,QC3name)
  c = plot_ethnic_roc_facut(s3_qc1, s3_qc2, s3_qc3, col,'Wightman',QC1name, QC2name,QC3name)
  prow <- plot_grid(a+ theme(legend.position="none"), 
                    b+ theme(legend.position="none",axis.text.y = element_blank())+ ylab(NULL), 
                    c+ theme(legend.position="none",axis.text.y = element_blank())+ylab(NULL),
                    ncol = 3, nrow = 1)
  legend_b <- get_legend(
    c+ guides(color = guide_legend(nrow = 1)) +
      theme(legend.position = "bottom")
  )
  plot_grid(prow, legend_b, ncol=1,rel_heights=c(3,.4))
  
}

#plot_auc_facut_all_sumstat(kunkle_adsp,kunkle_qc_variant_sumstat,kunkle_adsp_qc, bellenguez_adsp,  bellenguez_qc_on_variant_sumstat,bellenguez_adsp_qc, wightman_adsp,wightman_qc_variant_sumstat, wightman_adsp_qc)

plot_ethnic_roc_facut <- function(QC1, QC2, QC3, col, title,QC1name="QC_on_base", QC2name="QC_on_base+variants",QC3name="no_qc"){
  QC1 = plot_ethnic_roc(QC1, '',col, plot="F")
  QC1$qc_status = QC1name
  QC2 = plot_ethnic_roc(QC2, '', col, plot="F")
  QC2$qc_status = QC2name
  QC3 = plot_ethnic_roc(QC3, '', col, plot="F")
  QC3$qc_status = QC3name
  all = rbind(QC1,QC2, QC3)
  all %>% 
    filter(PRS !="PRS-001" & PRS != "PRS_1") 
  
  all$PRS = str_replace(all$PRS, 'PRS_', 'p = 0.')
  all[all$PRS=='p = 0.e5',]$PRS =" p = 1e-5"
  
  
  plot <- ggplot(data = all, aes(x=auc, y = PRS, shape=qc_status, color = qc_status))+
    geom_point(size=3)+
    facet_wrap(~ethnicity, ncol=1)+
    xlab('AUC')+ ggtitle(title)+
    theme_bw()
  return(plot)
}

plot_auc_facut_all_sumstat <- function(s1_qc1, s1_qc2, s1_qc3,s2_qc1, s2_qc2, s2_qc3,s3_qc1, s3_qc2, s3_qc3,QC1name="QC_on_base", QC2name="QC_on_base+variants",QC3name="no_qc", col = col_roc_E5){
  a = plot_ethnic_roc_facut(s1_qc1, s1_qc2, s1_qc3, col,'Kunkle',QC1name, QC2name,QC3name)
  b = plot_ethnic_roc_facut(s2_qc1, s2_qc2, s2_qc3, col,'Bellenguez',QC1name, QC2name,QC3name)
  c = plot_ethnic_roc_facut(s3_qc1, s3_qc2, s3_qc3, col,'Wightman',QC1name, QC2name,QC3name)
  prow <- plot_grid(a+ theme(legend.position="none"), 
                    b+ theme(legend.position="none",axis.text.y = element_blank())+ ylab(NULL), 
                    c+ theme(legend.position="none",axis.text.y = element_blank())+ylab(NULL),
                    ncol = 3, nrow = 1)
  legend_b <- get_legend(
    c+ guides(color = guide_legend(nrow = 1)) +
      theme(legend.position = "bottom")
  )
  plot_grid(prow, legend_b, ncol=1,rel_heights=c(3,.4))
  
}

## this will take awhile. ------
plot_R2_facut_allsumstat(kunkle_qc_base, kunkle_qc_variant_sumstat, kunkle_cT,bellenguez_qc_on_base, bellenguez_qc_on_variant_sumstat, bellenguez_cT0224, 
                         wightman_qc_base, wightman_qc_variant_sumstat, wightman_cT)

plot_auc_facut_all_sumstat(kunkle_qc_base, kunkle_qc_variant_sumstat, kunkle_cT,bellenguez_qc_on_base, bellenguez_qc_on_variant_sumstat,bellenguez_cT0224, 
                      wightman_qc_base, wightman_qc_variant_sumstat, wightman_cT) 

#plot_auc_facut_all_sumstat(kunkle_adsp,kunkle_qc_variant_sumstat,kunkle_adsp_qc, bellenguez_adsp,  bellenguez_qc_on_variant_sumstat,bellenguez_adsp_qc, wightman_adsp,wightman_qc_variant_sumstat, wightman_adsp_qc)
plot_R2_facut_allsumstat(kunkle_adsp,kunkle_qc_variant_sumstat,kunkle_adsp_qc, bellenguez_adsp,  bellenguez_qc_on_variant_sumstat,bellenguez_adsp_qc, wightman_adsp,wightman_qc_variant_sumstat, wightman_adsp_qc,
                         QC1name="No QC", QC2name="QC on variant and sumstat", QC3name="QC on all", boot_num=50)

##UKBB
plot_R2_facut_allsumstat(kunkle_adsp,kunkle_adsp_qc,kunkle_UKBB_qc, bellenguez_adsp,bellenguez_adsp_qc, bellenguez_UKBB_qc, wightman_adsp, wightman_adsp_qc,wightman_UKBB_qc,
                         QC1name="No QC",  QC2name="QC on all", QC3name="QC, UKBB variants only",boot_num=50)






## R2
## polyfun ------
mod1 <- glm(Diagnosis ~ Sex + Age + PRS, data= extract_afr(bellenguez_polypred_new), family=binomial)
RsqGLM(mod1, plot=FALSE)$Nagelkerke

#RsqGLM(glm(Diagnosis~PRS1+Sex+Age,family = 'binomial', data = bellenguez_polypred_new), plot=FALSE)$Nagelkerke
RsqGLM(glm(Diagnosis~PRS+Sex+Age,family = 'binomial', data = bellenguez_polypred_new), plot=FALSE)$Nagelkerke


plot_ethnic_R2(kunkle_polypred, "kunkle_polypred",  polypred_col, 50, "max_snp_per LD = ")
plot_ethnic_R2(wightman_polypred, "wightman_polypred",  polypred_col, 50, "max_snp_per LD = ")
plot_ethnic_R2(bellenguez_fixed_0224, "bellenguez_polypred",  polypred_col, 50, "max_snp_per LD = ")

## bellenguez------
#df,prs,main_title, plot = TRUE, legend="PRS_", boot_num = FALSE, replace="pT = 0."
plot_ethnic_R2(bellenguez_updateRSID, title, col, boot_num, replace)
plot_R2_cT(bellenguez_cT, bellenguez_qc, bellenguez_qc_on_base, bellenguez_qc_on_target, 'bellenguez_updateRSID', col_roc, plot=TRUE, boot_num = 50)
plot_R2_cT(bellenguez_cT0224, bellenguez_qc0224, bellenguez_qc_on_base0224, bellenguez_qc_on_target0224, 'bellenguez_fixed0224', col_roc_E5, boot_num=50 )

## plot beta and pvalue
p_beta_plot(bell_cT,c(4,8,12),'bellenguez_c+pT')
p_beta_plot(bell_log,c(4,8,12,16),'bellenguez_updated_polyfunPred')
p_beta_plot(bell_log_eur,c(4,8,12,16),'bellenguez_updated (EUR only)')
p_beta_plot(bell_susie_log,c(4,8,12,16),'bellenguez_susie')
p_beta_plot(bell_qc_log,c(4,8,12),'bellenguez_c+pT, QC')
p_beta_plot(bell_qc_EUR_log,c(4,8,12),'bellenguez_c+pT, QC, EUR')

## interested SNP
bellenguez_qc_interested
mod1 <- glm(Diagnosis ~ Sex + Age + PRS, data= bellenguez_interested, family=binomial)
 RsqGLM(mod1, plot=FALSE)$Nagelkerke

## kunkle----
## remove all APOE
plot_ethnic_R2(kunkle_withoutAPOE, col_roc_E5,"Kunkle without APOE", 5)
plot_ethnic_R2(kunkle_withoutAPOE_qc, col_roc_E5,"Kunkle without APOE, QC on sumstat_variant", 50)


RsqGLM(glm(Diagnosis~PRS1+Sex+Age,family = 'binomial', data = kunkle_polypred), plot=FALSE)$Nagelkerke
RsqGLM(glm(Diagnosis~PRS10+Sex+Age,family = 'binomial', data = kunkle_polypred), plot=FALSE)$Nagelkerke
 
## remove two APOE allele
col_roc_APOE = list("only2SNP","no2SNP","no2SNP_qc")
kunkle_cT_log<-log_reg(kunkle_APOE, col_roc_APOE, "kunkle_cT", plot=FALSE)   
kunkle_cT_qc_log_eur<-log_reg(extract_eur(kunkle_APOE), col_roc_APOE, "kunkle QC", plot=FALSE)
kunkle_cT_qc_log_afr<-log_reg(extract_afr(kunkle_APOE), col_roc_APOE, "kunkle QC", plot=FALSE)
kunkle_cT_qc_log_amr<-log_reg(kunkle_APOE[kunkle_APOE$final_population == "AMR",], col_roc_APOE, "kunkle QC", plot=FALSE)

##all population
plot_ethnic_R2(kunkle_cT, col_roc_E5, 'Kunkle', 50)
plot_ethnic_R2(kunkle_qc, col_roc_E5, 'Kunkle_QC, MAF 1%', 50)
plot_ethnic_R2(kunkle_qc_maf, col_roc_E5, 'Kunkle_QC, MAF 0.1%', 50)


##plot beta and p value
p_beta_plot(kunkle_cT_log,c(4,8,12),'kunkle')
p_beta_plot(kunkle_cT_qc_log,c(4,8,12),'kunkle_c+pT, QC')
p_beta_plot(kunkle_cT_qc_EUR_log,c(4,8,12),'kunkle_c+pT, QC, EUR')


## wightman ------
plot_R2_cT(wightman_cT, wightman_qc,wightman_qc_base,wightman_qc_target, 'Wightman', col_roc_E5, boot_num = 50)
plot_R2_cT(wightman_cT, wightman_qc,wightman_qc_base,wightman_qc_target, 'Wightman, EUR', col_roc_E5, boot_num = 50, eur=TRUE)
plot_ethnic_R2(wightman_cT, col_roc_E5, 'Wightman', 50)


##susie ------
plot_ethnic_R2(susie, col_roc_polypred, "bellenguez_susie",50)
plot_ethnic_R2(kunkle_susie, col_roc_polypred, "kunkle_susie",50)

RsqGLM(glm(Diagnosis~PRS10+Sex+Age,family = 'binomial', data = kunkle_susie), plot=FALSE)$Nagelkerke ## 0.14%
RsqGLM(glm(Diagnosis~PRS1+Sex+Age,family = 'binomial', data = kunkle_susie), plot=FALSE)$Nagelkerke ##0.139


##others -----
sbayesr_log <-  glm(Diagnosis~PRS+Sex+Age,family = 'binomial', data = sbayesR)
prsice_log <- glm(Diagnosis~PRS+Sex+Age,family = 'binomial', data = PRSice)


# Density
PRS_density <- function(df, name, Pthreshold){
  if(Pthreshold == 'PRS_e5'){
    pvalue=str_replace(Pthreshold,"PRS_e5",', p = 1e-5')
  }
  else{
    pvalue=str_replace(Pthreshold,"PRS_",', p = 0.')
  }

  ylimit = max(max(density(df[df$Diagnosis == 1,Pthreshold])$y), max(density(df[df$Diagnosis == 0,Pthreshold])$y))
  plot(density(df[df$Diagnosis == 0,Pthreshold]), main = paste(name,pvalue), col = "black", xlab = "PRS", ylim = c(0,ylimit*1.02))
  lines(density(df[df$Diagnosis == 1,Pthreshold]), col = "red")
}

plot_density <- function(df, df_qc, df_qc_base, df_qc_target, sumstat_name, Pthreshold , eur=FALSE){
  if(eur==TRUE){
    print("plot EUR subset")
    df = extract_eur(df)
    df_qc = extract_eur(df_qc)
    df_qc_base = extract_eur(df_qc_base)
    df_qc_target = extract_eur(df_qc_target)
  }
  PRS_density(df,  paste(sumstat_name, ", c+pT"), Pthreshold)
  PRS_density(df_qc,  paste(sumstat_name, ", QC"), Pthreshold)
  PRS_density(df_qc_base,  paste(sumstat_name, ", QC on summary stat"),Pthreshold)
  PRS_density(df_qc_target, paste(sumstat_name, ", QC on target"),Pthreshold)
}


PRS_density(wightman_qc_variant_sumstat, "Wightman, QC on variant_sumstat","PRS_5")
PRS_density(extract_eur(wightman_qc_variant_sumstat), "Wightman, EUR, QC on variant_sumstat","PRS_05")

plot_density(kunkle_cT, kunkle_qc, kunkle_qc_base, kunkle_qc_target, "kunkle", "PRS_5", eur=FALSE)
plot_density(kunkle_cT, kunkle_qc, kunkle_qc_base, kunkle_qc_target, "kunkle_eur", "PRS_05", eur=TRUE)
legend(-0.32,69, inset=.02, title="AD diagnosis",
       c("Case","Control"), fill=c("red","black"), horiz=TRUE, cex=0.8,xpd="NA")

plot_density(kunkle_cT, kunkle_qc, kunkle_qc_base, kunkle_qc_target, "kunkle_eur", "PRS_5", eur=TRUE)
legend(-0.12,179, inset=.02, title="AD diagnosis",
       c("Case","Control"), fill=c("red","black"), horiz=TRUE, cex=0.8,xpd="NA")

par(mfrow=c(2,2),xpd=FALSE)

PRS_density(extract_eur(kunkle_cT),"Kunkle, EUR","PRS_5")
PRS_density(extract_eur(kunkle_qc),"Kunkle, EUR, QC both","PRS_5")
#PRS_density(extract_eur(kunkle_qc_target_maf),"qc on target, MAF=0.1%,","PRS_5")
PRS_density(extract_eur(kunkle_qc_variant_sumstat), "Kunkle, EUR, QC on variant_sumstat","PRS_5")
PRS_density(extract_eur(kunkle_qc_variant), "QC on variants only, EUR,  MAF=0.1%", "PRS_5")
PRS_density(extract_eur(kunkle_qc_individual), "QC on individual only, EUR, MAF=0.1%", "PRS_5")
legend(-0.2,179, inset=.02, title="AD diagnosis",
       c("Case","Control"), fill=c("red","black"), horiz=TRUE, cex=0.8,xpd="NA")


plot_density(bellenguez_cT, bellenguez_qc, bellenguez_qc_on_base, bellenguez_qc_on_target, "bellenguez, EUR", "PRS_001", eur=TRUE)
PRS_density(extract_eur(bellenguez_cT), "bellenguez_EUR", "PRS_001")
PRS_density(extract_eur(bellenguez_qc_on_target), "EUR, qc on variant & individual","PRS_001")
PRS_density(extract_eur(bellenguez_qc_variant), "EUR, qc on variant", "PRS_001")
PRS_density(extract_eur(bellenguez_qc_individual), "EUR, qc on individual","PRS_001")

plot_density(wightman_cT, wightman_qc,wightman_qc_base,wightman_qc_target, "wightman", "PRS_e5")
legend(-2,1.1, inset=.02, title="AD diagnosis",
       c("Case","Control"), fill=c("red","black"), horiz=TRUE, cex=0.8,xpd="NA")

plot_density(wightman_cT, wightman_qc,wightman_qc_base,wightman_qc_target, "wightman", "PRS_e5", eur=TRUE)

PRS_density(extract_eur(wightman_qc), "wightman_EUR", "PRS_e5")
PRS_density(extract_eur(wightman_qc_target), "EUR, qc on variant & individual","PRS_e5")
PRS_density(extract_eur(wightman_qc_variant), "EUR, qc on variant", "PRS_e5")
PRS_density(extract_eur(wightman_qc_individual), "EUR, qc on individual","PRS_e5")


PRS_density(wightman_qc, "wightman, qc")
PRS_density(wightman_qc_target, "wightman, qc on target, p=0.5")
PRS_density(wightman_qc_variant, "wightman, qc on variant, p=0.5")
PRS_density(wightman_qc_individual, "wightman, qc on individual, p=0.5")

legend(-0.3,170, inset=.02, title="AD diagnosis",
       c("Case","Control"), fill=c("red","black"), horiz=TRUE, cex=0.8,xpd="NA")




