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
rename_preprocess <-function(df){
  df = read.delim(df)  
  names(df)[names(df) == 'age_covariate'] <- "Age"
  names(df)[names(df) == 'AD_status_final'] <- "Diagnosis"
  df = pre_process(df, FILE = TRUE)
  return (df)
}
pre_process <- function(df, FILE=FALSE){
  if(FILE==FALSE){
    df<- read.table(df,header=T,fill = T)
  }
  df$Age <- as.numeric(as.character(df$Age))
  print(paste("original=",  dim(df)[1], "rows"))
  df <- df %>%
    filter(Diagnosis != -1 & Age >= 65)
  print(paste("filtered=",  dim(df)[1], "rows"))
  return(df)
}
sbayesR = rename_preprocess("/gpfs/commons/home/tlin/output/prs/sbayesR.tsv")
  
bellenguez_update <- rename_preprocess('/gpfs/commons/home/tlin/output/prs/bellenguez_bellenguez_updateRSID_prs_PC.tsv')
bellenguez_fixed_0224 <- pre_process("/gpfs/commons/home/tlin/output/prs/bellenguez/fixed_0224/bellenguez_fixed_0224.tsv")
bellenguez_cT <- pre_process('/gpfs/commons/home/tlin/output/prs/bellenguez/bellenguez_pT_PRS_withPC.tsv')

bellenguez_fixed_0224 <-read.table("/gpfs/commons/home/tlin/output/prs/bellenguez/fixed_0224/bellenguez_fixed_0224.tsv", header = TRUE)

bellenguez_qc <- pre_process('/gpfs/commons/home/tlin/output/prs/bellenguez/bellenguez_qc_pT_PRS.tsv')
bellenguez_qc_on_target <- pre_process('/gpfs/commons/home/tlin/output/prs/bellenguez/qc_on_target.tsv')
bellenguez_qc_on_base <- pre_process('/gpfs/commons/home/tlin/output/prs/bellenguez/qc_on_base.tsv')

bellenguez_cT0224 <- pre_process('/gpfs/commons/home/tlin/output/prs/bellenguez/fixed_0224/bellenguez_no_qc.tsv')
bellenguez_qc0224 <- pre_process('/gpfs/commons/home/tlin/output/prs/bellenguez/fixed_0224/bellenguez_qc.tsv')
bellenguez_qc_on_target0224 <- pre_process('/gpfs/commons/home/tlin/output/prs/bellenguez/fixed_0224/bellenguez_qc_on_target.tsv')
bellenguez_qc_on_base0224 <- pre_process('/gpfs/commons/home/tlin/output/prs/bellenguez/fixed_0224/bellenguez_qc_on_base.tsv')

bellenguez_interested <- pre_process('/gpfs/commons/home/tlin/output/prs/bellenguez/updateRSID/interested_SNP/merged_updateRSID_interested_SNP.tsv') 
bellenguez_qc_interested <- pre_process('/gpfs/commons/home/tlin/output/prs/bellenguez/updateRSID/interested_SNP/merged_updateRSID_qc_interested_SNP.tsv')

##susie
susie <- rename_preprocess("/gpfs/commons/home/tlin/output/prs/bellenguez_susie/bellenguez_susie_prs_PC.tsv")

##kunkle
#kunkle_cT <- pre_process('/gpfs/commons/home/tlin/output/prs/kunkle/kunkle_pT_prs_before_qc.tsv')
kunkle_APOE <- pre_process("/gpfs/commons/home/tlin/output/prs/kunkle/fixed_0224/APOE.tsv")

##kunkle_qc
kunkle_cT <- pre_process('/gpfs/commons/home/tlin/output/prs/kunkle/fixed_0224/kunkle_no_qc.tsv')
kunkle_qc <- pre_process('/gpfs/commons/home/tlin/output/prs/kunkle/fixed_0224/kunkle_qc.tsv')
kunkle_qc_target <-pre_process('/gpfs/commons/home/tlin/output/prs/kunkle/fixed_0224/kunkle_qc_on_target.tsv')
kunkle_qc_base <- pre_process('/gpfs/commons/home/tlin/output/prs/kunkle/fixed_0224/kunkle_qc_on_base.tsv')

##PRSice
PRSice <- rename_preprocess("/gpfs/commons/home/tlin/output/prs/PRSice_pheno.tsv")

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

plot_density_diag <- function(df, col, title){
  for(i in 1:length(col)) {  
    plot(density((df[, col[i]])),col="grey", lty = 2, lwd = 2,main=str_replace(col[i],"PRS", "max_snp_per_LD = "), xlab='PRS')
    lines(density(df[df$Diagnosis == 1,col[i]]), col = 'red')
    lines(density(df[df$Diagnosis == 0, col[i]]), col = 'blue')
    legend("topright", legend=c("ALL","Case", "Control"),
           col=c("darkgrey","red", "blue"), lty=c(2,1,1), cex=0.8,
           box.lty=0)
  }
  mtext(title,                
        side = 3,line = -1.25,outer = TRUE)
}

par(mfrow=c(1,2))

plot_density_diag(bellenguez_update,c("PRS1","PRS10"), "bellenguez_update")
plot_density_diag(susie,c("PRS1","PRS10"), "susie")
bellenguez_eur = bellenguez_update[bellenguez_update$final_population == "EUR",]
plot_density_diag(bellenguez_eur,c("PRS1","PRS10"), "bellenguez_EUR")


par(mfrow=c(2,2))

plot_cs_density<- function(df, name){
  plot(density((df$PRS1)),col="red", lty = 1, lwd = 1,main=name, xlab='PRS', xlim = c(-0.004,0.003), ylim = c(0,1150))
  lines(density(df$PRS3), col='blue',lwd = 3)
  lines(density(df$PRS5), col='brown',lwd = 2,lty = 2)
  lines(density(df$PRS10), col='darkgreen', lty = 3, lwd = 1)
  legend("topright", legend=c("max SNP = 1","max SNP = 3","max SNP = 5","max SNP = 10"),
         col=c("red", "blue","brown","darkgreen"), lty=c(1,1,2,3), lwd = c(1,3,2,1), cex=0.8,
         box.lty=0)
}
plot_cs_density(susie, "susie")
plot_cs_density(bellenguez_all2, "polyfun")
lines(density(susie$PRS1), col='darkgreen', lty = 3, lwd = 1)



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





## Extract EUR only
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


##AUC ------
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

col_roc <- list("PRS_e5","PRS_001","PRS_005","PRS_01","PRS_05","PRS_1","PRS_5")

roc_result <-function(df, title=' ',column_for_roc){
  if('PRS_e5' %in% column_for_roc){
    roc_list <- roc(Diagnosis ~ PRS_e5+PRS_001+PRS_005+PRS_01+PRS_05+PRS_1+PRS_5, data = df)
  }else if ('PRS10' %in% column_for_roc){
    roc_list <- roc(Diagnosis ~ PRS1+PRS3+PRS5+PRS7+PRS10, data = df)
  }else{
    roc_list <- roc(Diagnosis ~ PRS_001+PRS_005+PRS_01+PRS_05+PRS_1+PRS_5, data = df)
  }
  plot_roc(roc_list,column_for_roc, title=title)
}



###c+T (before and after qc)-----
#### bellenguez-----
##bellenguez
col_roc <- list("PRS_001","PRS_005","PRS_01","PRS_05","PRS_1","PRS_5")
col_roc_E5 <- list("PRS_e5","PRS_001","PRS_005","PRS_01","PRS_05","PRS_1","PRS_5")
#bellenguez_cT_plot <- roc_result(bellenguez_cT, title = "bellenguez, clumping+pT", column_for_roc = col_roc_E5 )

##bellenguez_qc
bellenguez_ct_plot = roc_result(bellenguez_cT0224, title = "bellenguez, c+pT" , column_for_roc = col_roc_E5)
bellenguez_qc_plot = roc_result(bellenguez_qc0224, title = "bellenguez, c+pT, QC" , column_for_roc = col_roc_E5)
bellenguez_qc_on_base_plot = roc_result(bellenguez_qc_on_base0224, title = "bellenguez, c+pT, QC on summary_stats" , column_for_roc = col_roc_E5)
bellenguez_qc_on_target_plot = roc_result(bellenguez_qc_on_target0224, title = "bellenguez, c+pT, QC on target" , column_for_roc = col_roc_E5)

plot_grid(bellenguez_ct_plot, bellenguez_qc_plot,bellenguez_qc_on_base_plot,bellenguez_qc_on_target_plot,ncol = 2, nrow = 2)


## bellenguez_EUR
bellenguez_ct_plot = roc_result(extract_eur(bellenguez_cT0224), title = "bellenguez, c+pT, EUR" , column_for_roc = col_roc_E5)
bellenguez_qc_plot = roc_result(extract_eur(bellenguez_qc0224), title = "bellenguez, c+pT, QC,EUR" , column_for_roc = col_roc_E5)
bellenguez_qc_on_base_plot = roc_result(extract_eur(bellenguez_qc_on_base0224), title = "bellenguez, c+pT, QC on summary_stats,EUR" , column_for_roc = col_roc_E5)
bellenguez_qc_on_target_plot = roc_result(extract_eur(bellenguez_qc_on_target0224), title = "bellenguez, c+pT, QC on target,EUR" , column_for_roc = col_roc_E5)

plot_grid(bellenguez_ct_plot, bellenguez_qc_plot,bellenguez_qc_on_base_plot,bellenguez_qc_on_target_plot,ncol = 2, nrow = 2)


##bellenguez_qc_eur
bellenguez_qc_eur_plot = roc_result(extract_eur(bellenguez_qc0224),column_for_roc = col_roc_E5, title = "bellenguez, c+pT, QC, EUR" )
bellenguez_qc_afr_plot = roc_result(extract_afr(bellenguez_qc0224),column_for_roc = col_roc_E5, title = "bellenguez, c+pT,QC,AFR")
bellenguez_qc_amr_plot = roc_result(bellenguez_qc0224[bellenguez_qc0224$final_population == "AMR",],column_for_roc = col_roc_E5, title = "bellenguez, c+pT,QC,AMR")

plot_grid(bellenguez_qc_plot,bellenguez_qc_eur_plot,bellenguez_qc_afr_plot, bellenguez_qc_amr_plot,ncol =2, nrow = 2)

##bellenguez_interested_SNP
auc(roc(Diagnosis~PRS, data = bellenguez_interested))  # 0.5187
auc(roc(Diagnosis~PRS, data = extract_eur(bellenguez_interested))) ##0.5308
auc(roc(Diagnosis~PRS, data = extract_afr(bellenguez_interested))) #0.529
auc(roc(Diagnosis~PRS, data = extract_amr(bellenguez_interested))) #0.4991

auc(roc(Diagnosis~PRS, data = bellenguez_qc_interested)) # 0.5255
auc(roc(Diagnosis~PRS, data = extract_eur(bellenguez_qc_interested))) ##0.5422
auc(roc(Diagnosis~PRS, data = extract_afr(bellenguez_qc_interested))) #0.5215
auc(roc(Diagnosis~PRS, data = extract_amr(bellenguez_qc_interested))) #0.5182


####kunkle------

## Kunkle

col_roc_E5 <- list("PRS_e5","PRS_001","PRS_005","PRS_01","PRS_05","PRS_1","PRS_5")
### ori
kunkle_pt_plot = roc_result(extract_eur(kunkle_cT),"kunkle, clumping+pT, EUR", column_for_roc = col_roc_E5)

### Qced
kunkle_pt_qc_plot_update=roc_result(extract_eur(kunkle_qc),title="qc, clumping+pT, EUR", column_for_roc = col_roc_E5)

### base
kunkle_pt_qc_plot_base=roc_result(extract_eur(kunkle_qc_base), title="kunkle_qc_summary_stats,EUR",column_for_roc = col_roc_E5)

### target
kunkle_pt_qc_plot_target=roc_result(extract_eur(kunkle_qc_target),title="kunkle_qc_target,EUR", col_roc_E5)
plot_grid(kunkle_pt_plot,kunkle_pt_qc_plot_update,kunkle_pt_qc_plot_base,kunkle_pt_qc_plot_target,ncol = 2, nrow = 2)

## Kunkle_EUR
### ori
kunkle_pt_plot = roc_result(kunkle_cT,"kunkle, clumping+pT", column_for_roc = col_roc_E5)
kunkle_pt_plot_eur=roc_result( extract_eur(kunkle_cT),title="kunkle(EUR)",column_for_roc = col_roc_E5)
kunkle_pt_plot_afr=roc_result( extract_afr(kunkle_cT),title="kunkle(AFR)",column_for_roc = col_roc_E5)
kunkle_pt_amr_plot = roc_result(kunkle_cT[kunkle_cT$final_population == "AMR",], title = "kunkle(AMR)",column_for_roc = col_roc_E5)
plot_grid(kunkle_pt_plot,kunkle_pt_plot_eur,kunkle_pt_plot_afr, kunkle_pt_amr_plot,ncol = 2, nrow = 2)

### QCed EUR only
kunkle_pt_qc_plot_update=roc_result(extract_eur(kunkle_qc),title="qc, clumping+pT (EUR)",column_for_roc = col_roc_E5)

### base 
kunkle_pt_qc_plot_base=roc_result(extract_eur(kunkle_qc_base), title="kunkle_qc_summary_stats(EUR)",column_for_roc = col_roc_E5)

### target
kunkle_pt_qc_plot_target=roc_result(extract_eur(kunkle_qc_target),title="kunkle_qc_target(EUR)",column_for_roc = col_roc_E5)
plot_grid(kunkle_pt_plot,kunkle_pt_qc_plot_update,kunkle_pt_qc_plot_base,kunkle_pt_qc_plot_target,ncol = 2, nrow = 2)


## APOE
col_roc <- list("only2SNP","no2SNP","no2SNP_qc")
kunkle_APOE_roc <- roc(Diagnosis ~ only2SNP+no2SNP+no2SNP_qc, data = kunkle_APOE)
auc(roc(Diagnosis~only2SNP, data = kunkle_APOE)) ##0.621
auc(roc(Diagnosis~no2SNP, data = kunkle_APOE))##0.5131

auc(roc(Diagnosis~only2SNP, data = extract_eur(kunkle_APOE))) ##0.621
auc(roc(Diagnosis~no2SNP, data = extract_eur(kunkle_APOE)))  ##0.5834


### polyfun----

##bellenguez_all

bellenguez_roc = roc_result(bellenguez_fixed_0224, title="bellenguez(polyfun_pred)", column_for_roc = polypred_col)
bellenguez_eur_roc = roc_result(extract_eur(bellenguez_fixed_0224), title="bellenguez(polyfun_pred, EUR)", column_for_roc = polypred_col)
plot_grid(bellenguez_roc, bellenguez_eur_roc,ncol = 2, nrow = 1)

 ##susie
roc_list <- roc(Diagnosis ~ PRS1+PRS3+PRS5+PRS10, data = susie)
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
     



## Regression ----

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
        if(dim(R2)[1]==length(prs) && plot == TRUE){
          plotR2_boot(R2,main_title, prs,boot_num)
          return (R2)
        }
      }
    }  ##bootstrap
    else{
      mod <- glm(formula = frm,family = 'binomial', data = df)
      mod_R2 <-  RsqGLM(mod, plot=FALSE)$Nagelkerke 
      ## table
      if(i ==1){
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
        if(dim(cof)[1]/4==length(prs) && plot == TRUE){
          plotR2_boot(cof,main_title, prs,boot_num)
          return (cof)
        }
      }
    }
  }
}

log_reg(kunkle_cT, col_roc,"test", boot_num = 6)




plotR2 <- function(log_output, header){
  prs_col = seq(from=1, to=dim(log_output)[1], by=4)
  plot(log_output$R2_prs[prs_col], ylim  = c(min(log_output$R2[1], min(log_output$R2_prs[prs_col])), max(log_output$R2[1], max(log_output$R2_prs[prs_col]))),
       ylab='R-squared', xlab = "different PRS",pch=4,col="darkgreen", main=header,xaxt = "n")
  text(seq(from=1,to=length(prs_col)),log_output$R2_prs[prs_col] , 
       labels = log_output$R2_prs[prs_col], adj = c(0.5,2), xpd = TRUE, cex=1, col="darkgreen") 
  
  abline(h=log_output$NULL_R2[1], col=c("blue"), lty=2, lwd=4)
  text((length(prs_col)+1)/2,log_output$R2[1] , labels =  paste("Null model: ",log_output$R2[1]), adj = c(0.5,-1), xpd = TRUE, cex=1, col="blue")
  axis(1,                         # Define x-axis manually
       at = 1:length(prs_col),
       labels = c(str_replace(log_output$PRS[1],"PRS_e5","p = 1*e-5"),str_replace(log_output$PRS[seq(from=1, to=dim(log_output)[1], by=4)][-1],"PRS_","p = 0.")))
  
}

par(xpd = FALSE,mfrow=c(1,1))
plotR2_boot <- function(log_output, header, prs_col, boot){
  par(xpd = FALSE,mfrow=c(1,1))
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
plotR2_boot(test, "test", col_roc, FALSE)
#test_boot = log_reg(kunkle_cT, col_roc,"test", boot_num = 2)

#test = log_reg(kunkle_cT, col_roc,"test", plot = T)



par(mfrow=c(1,1)) 

rsq <- function(formula, data, indices) {
  d <- data[indices,] # allows boot to select sample
  fit <- glm(formula = formula,family = 'binomial', data = d)
  return(RsqGLM(fit, plot = FALSE)$Nagelkerke)
}
# bootstrapping with 1000 replications
results <- boot(data=kunkle_cT, statistic=rsq,
                R=20, formula=Diagnosis~Sex+Age,)

# view results
results
plot(results)

# get 95% confidence interval
boot.ci(results, type="bca")
result = boot.ci(boot.out=results,type='norm')
print(c(mean(results$t) - 2*sd(results$t), mean(results$t) + 2*sd(results$t)))




##bellenguez------
bell_log <- log_reg(bellenguez_update, col_roc_E5, 'bellenguez_updated', plot = FALSE, legend="PRS", replace="max_snp_per_locus=")
bell_log_eur <- log_reg(bellenguez_eur, col_roc_E5, 'bellenguez_updated', plot=FALSE, legend="PRS", replace="max_snp_per_locus=")
bell_susie_log<- log_reg(susie, list("PRS1","PRS3","PRS5",'PRS10'), 'bellenguez_susie', 'all',plot = FALSE, legend="PRS", replace="max_snp_per_locus=")

bell_cT_log<-log_reg(bellenguez_cT0224, col_roc_E5, "bellenguez_cT", plot=FALSE)   
bell_cT_qc_log<-log_reg(bellenguez_qc0224, col_roc_E5, "bellenguez QC", plot=FALSE)
bell_cT_qc_base_log<-log_reg(bellenguez_qc_on_base0224, col_roc, "bellenguez QC","qc on summary stats", plot=FALSE)
bell_cT_qc_target_log<-log_reg(bellenguez_qc_on_target0224, col_roc, "bellenguez QC on target", plot=FALSE)


par(mfrow=c(2,2),xpd=FALSE)
plotR2(bell_cT_log,"bellenguez")
plotR2(bell_cT_qc_log,"bellenguez, QC")
plotR2(bell_cT_qc_base_log,"bellenguez, QC_summary_stats")
plotR2(bell_cT_qc_target_log,"bellenguez, QC_base")


bell_cT_log<-log_reg(extract_eur(bellenguez_cT0224), col_roc_E5, "bellenguez_cT, EUR", plot=FALSE)   
bell_cT_qc_log<-log_reg(extract_eur(bellenguez_qc0224), col_roc_E5, "bellenguez QC, EUR", plot=FALSE)
bell_cT_qc_base_log<-log_reg(extract_eur(bellenguez_qc_on_base0224), col_roc, "bellenguez_QC_on summary stats, EUR", plot=FALSE)
bell_cT_qc_target_log<-log_reg(extract_eur(bellenguez_qc_on_target0224), col_roc, "bellenguez QC on target, EUR", plot=FALSE)


par(mfrow=c(2,2),xpd=FALSE)
plotR2(bell_cT_log,"bellenguez_cT, EUR")
plotR2(bell_cT_qc_log,"bellenguez, QC, EUR")
plotR2(bell_cT_qc_base_log,"bellenguez, QC_summary_stats, EUR")
plotR2(bell_cT_qc_target_log,"bellenguez, QC_base, EUR")


##kunkle----
## APOE
col_roc_APOE = list("only2SNP","no2SNP","no2SNP_qc")
kunkle_cT_log<-log_reg(kunkle_APOE, col_roc_APOE, "kunkle_cT", plot=FALSE)   
kunkle_cT_qc_log_eur<-log_reg(extract_eur(kunkle_APOE), col_roc_APOE, "kunkle QC", plot=FALSE)
kunkle_cT_qc_log_afr<-log_reg(extract_afr(kunkle_APOE), col_roc_APOE, "kunkle QC", plot=FALSE)
kunkle_cT_qc_log_amr<-log_reg(kunkle_APOE[kunkle_APOE$final_population == "AMR",], col_roc_APOE, "kunkle QC", plot=FALSE)


par(mfrow=c(2,2),xpd=FALSE)
plotR2(kunkle_cT_log,"kunkle, APOE")
plotR2(kunkle_cT_qc_log_eur,"kunkle, EUR")
plotR2(kunkle_cT_qc_log_afr,"kunkle, AFR")
plotR2(kunkle_cT_qc_log_amr,"kunkle, AMR")



##all population
kunkle_cT_log<-log_reg(kunkle_cT, col_roc_E5, "kunkle_cT", plot=FALSE)   
kunkle_cT_qc_log<-log_reg(kunkle_qc, col_roc_E5, "kunkle QC", plot=FALSE)
kunkle_cT_qc_base_log<-log_reg(kunkle_qc_base, col_roc, "kunkle QC","qc on summary stats", plot=FALSE)
kunkle_cT_qc_target_log<-log_reg(kunkle_qc_target, col_roc, "kunkle QC on target", plot=FALSE)

par(mfrow=c(2,2),xpd=FALSE)
plotR2(kunkle_cT_log,"kunkle, cT")
plotR2(kunkle_cT_qc_log,"kunkle, qc")
plotR2(kunkle_cT_qc_base_log,"kunkle, qc on summary stats")
plotR2(kunkle_cT_qc_target_log,"kunkle, qc on target")

## EUR
kunkle_cT_log<-log_reg(extract_eur(kunkle_cT), col_roc_E5, "kunkle_cT, EUR", plot=FALSE)   
kunkle_cT_qc_log<-log_reg(extract_eur(kunkle_qc), col_roc_E5, "kunkle QC,EUR", plot=FALSE)
kunkle_cT_qc_base_log<-log_reg(extract_eur(kunkle_qc_base), col_roc, "kunkle QC","qx on summary stats", plot=FALSE)
kunkle_cT_qc_target_log<-log_reg(extract_eur(kunkle_qc_target), col_roc, "kunkle QC on target,EUR", plot=FALSE)


par(mfrow=c(2,2))
plotR2(kunkle_cT_log,"kunkle, cT, EUR")
plotR2(kunkle_cT_qc_log,"kunkle, qc, EUR")
plotR2(kunkle_cT_qc_base_log,"kunkle, qc on summary stats, EUR")
plotR2(kunkle_cT_qc_target_log,"kunkle, qc on target, EUR")

##other
kunkle_cT_log_amr <- log_reg(kunkle_cT[kunkle_cT$final_population == "AMR",], col_roc_E5, "kunkle_cT, AMR", plot=FALSE)
kunkle_cT_log_afr<-log_reg(extract_afr(kunkle_cT), col_roc_E5, "kunkle_cT, AFR", plot=FALSE)   

par(mfrow=c(1,2))
plotR2(kunkle_cT_log_afr,"kunkle, cT, AFR")
plotR2(kunkle_cT_log_amr,"kunkle, cT, AMR")
plotR2(kunkle_cT_log,"kunkle, cT, EUR")

sbayesr_log <-  glm(Diagnosis~PRS+Sex+Age,family = 'binomial', data = sbayesR)
prsice_log <- glm(Diagnosis~PRS+Sex+Age,family = 'binomial', data = PRSice)


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
p_beta_plot(bell_cT,c(4,8,12),'bellenguez_c+pT')
p_beta_plot(bell_log,c(4,8,12,16),'bellenguez_updated_polyfunPred')
p_beta_plot(bell_log_eur,c(4,8,12,16),'bellenguez_updated (EUR only)')
p_beta_plot(bell_susie_log,c(4,8,12,16),'bellenguez_susie')
p_beta_plot(kunkle_cT_log,c(4,8,12),'kunkle')



p_beta_plot(kunkle_cT_qc_log,c(4,8,12),'kunkle_c+pT, QC')
p_beta_plot(kunkle_cT_qc_EUR_log,c(4,8,12),'kunkle_c+pT, QC, EUR')

p_beta_plot(bell_qc_log,c(4,8,12),'bellenguez_c+pT, QC')
p_beta_plot(bell_qc_EUR_log,c(4,8,12),'bellenguez_c+pT, QC, EUR')




NagelkerkeR2(mod2)$R2
RsqGLM(mod2)$Nagelkerke
par(mfrow=c(2,2),xpd=FALSE)


kunkle_eur = extract_eur(kunkle_cT)
plot(density(kunkle_eur[kunkle_eur$Diagnosis == 1,]$PRS), main = "kunkle_EUR", col = "red", xlab = "PRS")
lines(density(kunkle_eur[kunkle_eur$Diagnosis == 0,]$PRS), col = "black")


extract_eur(kunkle_qc_base)
plot(density(extract_eur(kunkle_qc_base)[extract_eur(kunkle_qc_base)$Diagnosis == 1,]$PRS), main = "kunkle_EUR_QC", col = "red", xlab = "PRS")
lines(density(extract_eur(kunkle_qc_base)[extract_eur(kunkle_qc_base)$Diagnosis == 0,]$PRS), col = "black")

PRS_density <- function(df, name){
  plot(density(df[df$Diagnosis == 0,]$PRS_5), main = name, col = "black", xlab = "PRS")
  lines(density(df[df$Diagnosis == 1,]$PRS_5), col = "red")
}

PRS_density(extract_eur(kunkle_cT),"Kunkle, EUR (P = 0.5)")
PRS_density(extract_eur(kunkle_qc),"Kunkle, qc, EUR (P = 0.5)")
PRS_density(extract_eur(kunkle_qc_base),"qc on summary stats, EUR (P = 0.5)")
PRS_density(extract_eur(kunkle_qc_target),"qc on target, EUR (P = 0.5)")

legend(-0.25,80, inset=.02, title="AD diagnosis",
       c("Case","Control"), fill=c("red","black"), horiz=TRUE, cex=0.8,xpd="NA")
