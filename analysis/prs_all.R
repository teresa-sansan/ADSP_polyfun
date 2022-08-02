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
library("gridGraphics")

## load PRS data ----
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
### kunkle ----
kunkle_APOE <- pre_process("/gpfs/commons/home/tlin/output/prs/new_plink/kunkle/APOE_SNP_qc.tsv")
kunkle_withoutAPOE_qc <- pre_process("/gpfs/commons/home/tlin/output/prs/new_plink/kunkle/kunkle_ADSP_no_apoe.tsv")

kunkle_adsp <- pre_process('/gpfs/commons/home/tlin/output/prs/new_plink/kunkle/kunkle_ADSP.tsv')
kunkle_adsp_qc <- pre_process('/gpfs/commons/home/tlin/output/prs/new_plink/kunkle/kunkle_ADSP_qc_all.tsv')
kunkle_qc_maf <- pre_process('/gpfs/commons/home/tlin/output/prs/kunkle/fixed_0224/kunkle_qc_all_maf01.tsv')

kunkle_qc_target_maf <- pre_process('/gpfs/commons/home/tlin/output/prs/kunkle/fixed_0224/kunkle_qc_target_maf01.tsv')
kunkle_new_beta <- pre_process('/gpfs/commons/home/tlin/output/prs/kunkle/fixed_0224/new_beta_noqc.tsv')


### bellenguez ----- 
bellenguez_cT <- pre_process('/gpfs/commons/home/tlin/output/prs/bellenguez/fixed_0224/bellenguez_no_qc.tsv')
bellenguez_qc_all <- pre_process('/gpfs/commons/home/tlin/output/prs/new_plink/bellenguez/bellenguez_ADSP_qc_all.tsv')
bellenguez_qc_on_base <- pre_process('/gpfs/commons/home/tlin/output/prs/bellenguez/fixed_0224/bellenguez_qc_on_base.tsv')

bellenguez_interested <- pre_process('/gpfs/commons/home/tlin/output/prs/bellenguez/updateRSID/interested_SNP/merged_updateRSID_interested_SNP.tsv') 
bellenguez_qc_interested <- pre_process('/gpfs/commons/home/tlin/output/prs/bellenguez/updateRSID/interested_SNP/merged_updateRSID_qc_interested_SNP.tsv')
bellenguez_interest_max <- pre_process('/gpfs/commons/home/tlin/output/prs/bellenguez/updateRSID/interested_SNP/maxPIP/merged_updateRSID_qc_interested_SNP.tsv')
bellenguez_interest_min <- pre_process('/gpfs/commons/home/tlin/output/prs/bellenguez/updateRSID/interested_SNP/minPIP/merged_updateRSID_qc_interested_SNP.tsv')

bellenguez_qc_variant <- pre_process('/gpfs/commons/home/tlin/output/prs/bellenguez/fixed_0224/bellenguez_qc_on_variant.tsv')

### wightman ----
wightman_cT <- pre_process('/gpfs/commons/home/tlin/output/prs/wightman/before_qc.tsv')
wightman_qc <- pre_process('/gpfs/commons/home/tlin/output/prs/new_plink/wightman/wightman_ADSP_qc_all.tsv')
wightman_qc_base <- pre_process('/gpfs/commons/home/tlin/output/prs/wightman/qc_on_base.tsv')
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

bellenguez_adsp <- pre_process('/gpfs/commons/home/tlin/output/prs/new_plink/bellenguez/bellenguez_ADSP.tsv')
wightman_adsp <- pre_process('/gpfs/commons/home/tlin/output/prs/new_plink/wightman/wightman_ADSP.tsv')


bellenguez_adsp_qc <- pre_process('/gpfs/commons/home/tlin/output/prs/new_plink/bellenguez/bellenguez_ADSP_qc_all.tsv')
wightman_adsp_qc <- pre_process('/gpfs/commons/home/tlin/output/prs/new_plink/wightman/wightman_ADSP_qc_all.tsv')

kunkle_UKBB_qc <- pre_process('/gpfs/commons/home/tlin/output/prs/new_plink/kunkle/kunkle_ADSP_UKBB_qc.tsv')
bellenguez_UKBB_qc <- pre_process('/gpfs/commons/home/tlin/output/prs/new_plink/bellenguez/bellenguez_ADSP_UKBB_qc.tsv')
wightman_UKBB_qc <- pre_process('/gpfs/commons/home/tlin/output/prs/new_plink/wightman/wightman_ADSP_UKBB_qc.tsv')

### susie -----
kunkle_susie <- pre_process("/gpfs/commons/home/tlin/output/prs/polypred/kunkle/susie.prs.tsv")
susie <- pre_process("/gpfs/commons/home/tlin/output/prs/bellenguez_susie/bellenguez_susie_prs_PC.tsv")

## Others
PRSice <- pre_process("/gpfs/commons/home/tlin/output/prs/PRSice_pheno.tsv")
sbayesR = pre_process("/gpfs/commons/home/tlin/output/prs/sbayesR.tsv")

## Density plot
## plotting case/control plot ---- 
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

## create a race-separated pivot df for poly_pred (which is easier for plotting). 
## can be apply to different methods (polyfunpred or pT)
pivot_df<- function(df){ 
  df_new <- df %>%
    subset(final_population ==  "EUR" | final_population == "AMR"| final_population =="AFR") %>%   ## only keep the race of interest
    pivot_longer(                           
      cols = starts_with("PRS"),
      names_to = 'method', 
      values_to = 'PRS') %>% 
    as.data.frame()
  #df_new$Diagnosis = as.factor(df_new$Diagnosis)
  df_new$Diagnosis_state=ifelse(df_new$Diagnosis == 0, "control","case")  ## create a new column so that it will show case/control when plotting
  return(df_new)
}

## draw facet dist plot (by method) 
## can be apply to different methods (polyfunpred or pT)
##W here method = how many causal SNPs are assumed in the LD block/ pT threshold.
plot_facet_dist <- function(df, title){
  ggplot(pivot_df(df), aes(x =PRS, group=Diagnosis_state,fill = Diagnosis_state, color = Diagnosis_state))+
    geom_density(alpha = 0.01) +
    facet_grid(method~final_population)+theme_bw()+ggtitle(title)
} 

## plot race/method separated dist plots
plot_facet_dist(kunkle_polypred, "kunkle_polypred")
plot_faect_dist(kunkle_new_beta, "kunkle_new_beta")
plot_facet_dist(kunkle_susie, "kunkle_susie")
plot_facet_dist(wightman_beta, "wightman_new_beta")
plot_facet_dist(wightman_susie_max10_polypred, "wightman_susie_new_beta")


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
col_APOE <- list("PRS")
## AUC function ------
plot_roc <- function(roclist, roccol,title = FALSE){
  legendname = list()
  for (i in roccol){                
    auc_value = round(roclist[[i]]$auc,3)
    if(i == "PRS_e5"){
      legendname = append(legendname,paste('p = 1e-5, auc =',auc_value))
    }else if(!grepl("_",i, fixed = T)){
      legendname = append(legendname,paste(str_replace(i,"PRS",'max_snp_'),', auc =',auc_value))
    }else{
      legendname = append(legendname,paste(str_replace(i,"PRS_","p = 0."),', auc =',auc_value))
    }
  }
  curve <- ggroc(roclist)
  curve + ggtitle(title)+theme_bw()+
    scale_colour_discrete(name="PRS",labels = c(legendname))
}   

## A function that can calculate AUC with different method/threshold
## If plot = T(default), will return a plot with AUC using different method/threshold, if false, will return a list for 
## other purpose. 

roc_result <-function(df, title=' ',column_for_roc = col_roc_E5, plot = TRUE){
  if('PRS_e5' %in% column_for_roc){
    roc_list <- roc(Diagnosis ~ PRS_e5+PRS_001+PRS_005+PRS_01+PRS_05+PRS_1+PRS_5, data = df, quiet=T)
  }else if ('PRS3' %in% column_for_roc){
    roc_list <- roc(Diagnosis ~ PRS1+PRS3+PRS5+PRS7+PRS10, data = df, quiet=T)
  }else if ('PRS10' %in% column_for_roc){
    roc_list <- roc(Diagnosis ~ PRS1+PRS5+PRS10, data = df, quiet=T)
  }else if ('new_beta' %in% column_for_roc){
    roc_list <- roc(Diagnosis ~ new_beta + new_beta_cT +new_beta_susie, data=df, quiet=T)
  }else 
    roc_list <- roc(Diagnosis ~ PRS_001+PRS_005+PRS_01+PRS_05+PRS_1+PRS_5, data = df, quiet = T)
  
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
roc_result_boot <-function(df, title=' ',column_for_roc=col_roc_E5, boot_num=50, mean=FALSE, plot=FALSE){
  CI_roc = data.frame(matrix(ncol = 4, nrow = 0))
  colnames(CI_roc) = c("PRS","boot_CI_lower","boot_mean","boot_CI_upper")
  for (i in 1:length(column_for_roc)){
    col = column_for_roc[[i]]
    roc_formula <- roc(df[["Diagnosis"]],df[[col]], quiet=T)
    ci = ci.auc(roc_formula, conf.level = 0.95, method='bootstrap', boot.n = boot_num)
    CI_roc[nrow(CI_roc) + 1,] = c(col,ci)
  }
  ## only return auc if mean = TRUE
  if (mean==TRUE){ 
    roclist=list()
    for (i in 1:dim(CI_roc)[1]){
      roclist[i] = as.numeric(CI_roc$boot_mean[i])
    }
    return(roclist)
  }
  return(CI_roc)
}


# A function to plot result of diff race all at once. 
plot_ethnic_roc <- function(df, title, col,boot=TRUE,boot_num=5, plot=TRUE){
  if (boot != TRUE){ ## no bootstraping, use the OG roc_result
    EUR = roc_result(extract_eur(df), title = paste(title, ", EUR") , column_for_roc = col, plot=plot)
    AFR = roc_result(extract_afr(df), title = paste(title, ", AFR") , column_for_roc = col, plot=plot)
    AMR = roc_result(extract_amr(df), title = paste(title, ", AMR") , column_for_roc = col, plot=plot)
    if (plot == TRUE){
      plot_grid(EUR, AFR, AMR,ncol = 1, nrow = 3)
      return()
    }
  } else {
    EUR = roc_result_boot(extract_eur(df), column_for_roc = col, boot_num = boot_num)
    AFR = roc_result_boot(extract_afr(df), column_for_roc = col, boot_num = boot_num)
    AMR = roc_result_boot(extract_amr(df), column_for_roc = col, boot_num = boot_num)
  }
  df <- data_frame(matrix(ncol = 0, nrow = length(col)*3))
  df$PRS =  rep(unlist(col),3)
  df$ethnicity = c(rep("EUR",length(col)),rep("AFR",length(col)),rep("AMR",length(col)))
  df$ethnicity <- factor(df$ethnicity,      # Reordering group factor levels
                         levels = c("EUR","AFR","AMR"))
  if (boot != TRUE){
    df$auc = append(append(unlist(EUR),unlist(AFR)),unlist(AMR))
    return(df[,c(2:4)])
  }else{
    df$auc = append(append(unlist(EUR$boot_mean),unlist(AFR$boot_mean)),unlist(AMR$boot_mean))
    df$boot_CI_lower = append(append(unlist(EUR$boot_CI_lower),unlist(AFR$boot_CI_lower)),unlist(AMR$boot_CI_lower))
    df$boot_CI_upper = append(append(unlist(EUR$boot_CI_upper),unlist(AFR$boot_CI_upper)),unlist(AMR$boot_CI_upper))
    
    df[c("auc","boot_CI_lower","boot_CI_upper")] <- sapply(df[c("auc","boot_CI_lower","boot_CI_upper")],as.numeric)
    return(df[,c(2:6)])
    
  }
}
plot_ethnic_roc_facut <- function(QC1, QC2, QC3=FALSE, col=col_roc_E5, title=' ',QC1name="QC_on_base", QC2name="QC_on_base+variants",QC3name="no_qc",boot=TRUE, boot_num=50, legendname=FALSE){
  QC1 = plot_ethnic_roc(QC1, '', col, plot="F",boot=boot,boot_num=boot_num)
  QC1$qc_status = QC1name
  QC2 = plot_ethnic_roc(QC2, '', col, plot="F",boot=boot,boot_num=boot_num)
  QC2$qc_status = QC2name
  if(QC3 != FALSE){
    QC3 = plot_ethnic_roc(QC3, '', col, plot="F",boot=boot,boot_num=boot_num)
    QC3$qc_status = QC3name
    all = rbind(QC1,QC2, QC3)
  }else
    all  = rbind(QC1,QC2)
 
  all %>% 
    filter(PRS !="PRS-001" & PRS != "PRS_1") 
  all$PRS = str_replace(all$PRS, 'PRS_', 'p = 0.')
  all[all$PRS=='p = 0.e5',]$PRS =" p = 1e-5"
  plot <- ggplot(data = all, aes(x=auc, y = PRS, color = qc_status))+
    geom_point(size=3,alpha=0.8)+
    facet_wrap(~ethnicity, ncol=1)+
    xlab('AUC')+ ggtitle(title)+xlim(0.45, 0.75)+
    theme_bw()
  
  if (boot == TRUE){  ## plot error bar or not
    plot <- plot + geom_errorbar(aes(xmin=boot_CI_lower, xmax=boot_CI_upper), width=.2,alpha=0.5) 
  }
  if(legendname != FALSE){ 
    plot = plot +guides(col=guide_legend(legendname))
  }
  return(plot)
}



plot_auc_facut_all_sumstat <- function(s1_qc1, s1_qc2, s1_qc3,s2_qc1, s2_qc2, s2_qc3,s3_qc1, s3_qc2, s3_qc3, col = col_roc_E5, QC1name="QC_on_base", QC2name="QC_on_base+variants",QC3name="no_qc",boot=TRUE, boot_num=5){
  a = plot_ethnic_roc_facut(s1_qc1, s1_qc2, s1_qc3, col,'Kunkle',QC1name, QC2name,QC3name,boot=boot,boot_num=boot_num)
  b = plot_ethnic_roc_facut(s2_qc1, s2_qc2, s2_qc3, col,'Bellenguez',QC1name, QC2name,QC3name,boot=boot,boot_num=boot_num)
  c = plot_ethnic_roc_facut(s3_qc1, s3_qc2, s3_qc3, col,'Wightman',QC1name, QC2name,QC3name,boot=boot,boot_num=boot_num)
  
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


## every race should have diff APOE value
## fixed legend
plot_ethnic_roc_facut_add <- function(QC1, QC2, QC3, col=col_roc_E5, title=' ',QC1name="QC_on_base", QC2name="QC_on_base+variants",QC3name="no_qc",boot=TRUE, boot_num=50, legendname=FALSE, APOE=FALSE){
  
  QC1 = plot_ethnic_roc(QC1, '', col, plot="F",boot=boot,boot_num=boot_num)
  QC1$qc_status = QC1name
  QC2 = plot_ethnic_roc(QC2, '', col, plot="F",boot=boot,boot_num=boot_num)
  QC2$qc_status = QC2name
  if(class(QC3) != "logical"){
    QC3 = plot_ethnic_roc(QC3, '', col, plot="F",boot=boot,boot_num=boot_num)
    QC3$qc_status = QC3name
    all = rbind(QC1,QC2, QC3)
  }else
    all  = rbind(QC1,QC2)
  
  all %>% 
    filter(PRS !="PRS-001" & PRS != "PRS_1") 
  all$PRS = str_replace(all$PRS, 'PRS_', 'p = 0.')
  all[all$PRS=='p = 0.e5',]$PRS =" p = 1e-5"
  
  if (class(APOE) != "logical"){
    roc_formula= roc(APOE[["Diagnosis"]],APOE[["PRS"]],quiet=TRUE)
    if(boot==TRUE){
      APOE_CI = ci.auc(roc_formula, conf.level = 0.95, method='bootstrap', boot.n = boot_num )[1:3]
      #print(APOE_CI)
    } else{ ## no bootstrapping
      APOE_roc = roc_formula$auc
      print(roc_formula$auc)
    } 
  }
  
  plot <- ggplot(data = all, aes(x=auc, y = PRS, color = qc_status))+
    geom_point(size=3,alpha=0.8)+
    facet_wrap(~ethnicity, ncol=1)+
    xlab('AUC')+ ggtitle(title)+xlim(0.45, 0.75)+
    theme_bw()
 
  if (boot == TRUE){  ## plot error bar or not
    plot <- plot + geom_errorbar(aes(xmin=boot_CI_lower, xmax=boot_CI_upper), width=.2,alpha=0.5) 
  }
  if(legendname != FALSE){ 
    plot = plot +guides(col=guide_legend(legendname))
  }
  if(class(APOE) != "logical"){
    if(boot == TRUE){
      plot = plot + geom_vline(xintercept=APOE_CI[2],lwd=0.8,colour="black") + 
        geom_vline(xintercept=APOE_CI[1],lwd=0.5,colour="grey") +
        geom_vline(xintercept=APOE_CI[3],lwd=0.5,colour="grey")
    }else
      plot = plot + geom_vline(xintercept=APOE_roc,lwd=0.8,colour="black")
  }
  return(plot)
}



plot_ethnic_roc_facut_add(kunkle_qc_base, kunkle_qc_variant_sumstat, kunkle_cT,boot_num=2, APOE=kunkle_APOE)



## R2 functions ----
## calcualte pseudo rsquare 
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
  mod1 <- glm(Diagnosis ~ Sex + Age, data= df, family=binomial) ## first create a formula only using covariates
  mod1_R2 <- RsqGLM(mod1, plot=FALSE)$Nagelkerke
  for (i in 1:length(prs)){
    frm <- as.formula(paste("Diagnosis ~ Sex + Age+", prs[[i]])) ## add the PRS you wanted
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
    }  ##if don't do boostrap
    else{
      mod <- glm(formula = frm,family = 'binomial', data = df) ## with PRS
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
        if(dim(cof)[1]/4==length(prs) ){
          if (plot == TRUE)
            plotR2_boot(cof,main_title, prs,boot_num)
          return (cof)
        }
      }
    }
  }
}

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
  df = rbind(EUR, AFR, AMR)
  df$ethnicity <- factor(df$ethnicity,      # Reordering group factor levels
                         levels = c("EUR","AFR","AMR"))
  if(plot != "TRUE"){
    return(df)
  }
}

## plot facet 
## plot ethnic facut plot for one sumstat
plot_ethnic_R2_facut <- function(QC1, QC2, QC3=FALSE, col=col_roc_E5, boot_num=50, title=' ',QC1name="QC_on_base", QC2name="QC_on_base+variants",QC3name="no_qc", legendname=FALSE){
  QC1 = plot_ethnic_R2(QC1, col, '', boot_num, plot="F")
  QC2 = plot_ethnic_R2(QC2, col, '', boot_num, plot="F")
  QC1$qc_status = QC1name
  QC2$qc_status = QC2name
  if(QC3 != FALSE){
    QC3 = plot_ethnic_R2(QC3, col, '', boot_num, plot="F")
    QC3$qc_status = QC3name
    all = rbind(QC1,QC2, QC3)
  }else
    all = rbind(QC1, QC2)
  all %>% 
    filter(PRS !="PRS-001" & PRS != "PRS_1") 
  
  all$PRS = str_replace(all$PRS, 'PRS_', 'p = 0.')
  all[all$PRS=='p = 0.e5',]$PRS =" p = 1e-5"
  all$boot_mean = all$boot_mean * 100
  plot <- ggplot(data = all, aes(x= boot_mean, y = PRS, color = qc_status))+
    geom_point(size=3, alpha=0.75)+
    geom_pointrange(aes(xmin=boot_CI_lower*100, xmax=boot_CI_upper*100), linetype="dotted") +
    facet_wrap(~ethnicity, ncol=1)+
    xlab('R squared (%)')+ ggtitle(title)+xlim(-1, 25)+
    theme_bw()
  if(legendname != FALSE){  ## testout new legend title
    plot = plot +guides(col=guide_legend(legendname))
  }
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


## Regression (plot beta and pvalue)
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

##result ----
# c+T (before and after qc)-----
## bellenguez-----

## bellenguez_ethnics 
plot_ethnic_roc_facut(bellenguez_cT0224,bellenguez_adsp,bellenguez_adsp_qc,col_roc_E5,title='bellenguez',
                      QC1name= "k=old plink", QC2name="new_plink",QC3name='new_plink, qc')

##bellenguez_interested_SNP 
df_list <- list(bellenguez_interested, bellenguez_interest_max,bellenguez_interest_min,
                extract_eur(bellenguez_interested),extract_eur(bellenguez_interest_max),extract_eur(bellenguez_interest_min),
                extract_afr(bellenguez_interested),extract_afr(bellenguez_interest_max),extract_afr(bellenguez_interest_min),
                extract_amr(bellenguez_interested),extract_amr(bellenguez_interest_max),extract_amr(bellenguez_interest_min))

for (i in df_list){
  print(auc(roc(Diagnosis~PRS, data = i, quiet=T))[1])
}


## kunkle------

## Kunkle Qced

plot_ethnic_roc_facut(kunkle_qc_all, kunkle_withoutAPOE_qc,title='kunkle (QC)',
                      QC1name = "kunkle", QC2name = "kunkle_no_APOE", 
                      legendname = "APOE status", boot_num=2)


kunkle_pt_qc_plot_maf=roc_result(kunkle_qc_maf,title="QC, MAF >0 0.1%", column_for_roc = col_roc_E5)
##kunkle_pt_qc_plot_base=roc_result(kunkle_qc_base, title="kunkle_qc_summary_stats",column_for_roc = col_roc_E5)
plot_grid(kunkle_pt_qc_plot,kunkle_pt_qc_plot_maf,ncol = 2, nrow = 2)



##cros group 
par(mfrow=c(3,1))
col_roc_polypred <- list("PRS1","PRS3","PRS5","PRS7","PRS10")
plot_ethnic_roc (kunkle_cT, 'kunkle', col_roc_E5)
plot_ethnic_roc(kunkle_qc_variant_sumstat, "kunkle_qc_sumstat_variant", col_roc_E5)
plot_ethnic_roc (kunkle_polypred, 'kunkle', col_roc_polypred)
plot_ethnic_roc (bellenguez_fixed_0224,'bellenguez', col_roc_polypred)
plot_ethnic_roc(kunkle_new_beta, 'kunkle, (effectsize * pip) ',col_roc_E5)
plot_ethnic_roc(kunkle_susie, "kunkle(susie)", col_roc_polypred)
plot_ethnic_roc(susie, "bellenguez(susie)", col_roc_polypred)


## APOE
plot_ethnic_roc(kunkle_withoutAPOE, "kunkle remove APOE region", col_roc_E5)
col_roc <- list("only2SNP","no2SNP","no2SNP_qc")


kunkle_APOE_roc <- roc(Diagnosis ~ only2SNP+no2SNP+no2SNP_qc, data = kunkle_APOE)
auc(roc(Diagnosis~only2SNP, data = kunkle_APOE)) ##0.621
auc(roc(Diagnosis~only2SNP, data = extract_eur(kunkle_APOE))) ##0.6187

for (i in list(extract_eur(kunkle_APOE),extract_afr(kunkle_APOE),extract_amr(kunkle_APOE))){
  print(auc(roc(Diagnosis~no2SNP, data = i, quiet=T))[1])
}

plot_ethnic_roc_facut(kunkle_cT,kunkle_withoutAPOE,kunkle_withoutAPOE_qc,col_roc_E5, title='kunkle',
                      QC1name= "with APOE", QC2name= "without APOE", QC3name="without APOE, QC")

### wightman 
plot_ethnic_roc (wightman_polypred, 'wightman', col_roc_polypred)
plot_auc_cT(wightman_cT,wightman_qc, wightman_qc_base, wightman_qc_target, 'Wightman', col_roc_E5)
plot_ethnic_roc(wightman_cT, "wightman_qc", col_roc_E5)
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
     


## auc ------------


## this will take awhile. ------

plot_auc_facut_all_sumstat(kunkle_qc_base, kunkle_qc_variant_sumstat, kunkle_cT,bellenguez_qc_on_base, bellenguez_qc_on_variant_sumstat,bellenguez_cT0224, 
                      wightman_qc_base, wightman_qc_variant_sumstat, wightman_cT) 

#plot_auc_facut_all_sumstat(kunkle_adsp,kunkle_qc_variant_sumstat,kunkle_adsp_qc, bellenguez_adsp,  bellenguez_qc_on_variant_sumstat,bellenguez_adsp_qc, wightman_adsp,wightman_qc_variant_sumstat, wightman_adsp_qc)
plot_R2_facut_allsumstat(kunkle_adsp,kunkle_qc_variant_sumstat,kunkle_adsp_qc, bellenguez_adsp,  bellenguez_qc_on_variant_sumstat,bellenguez_adsp_qc, wightman_adsp,wightman_qc_variant_sumstat, wightman_adsp_qc,
                         QC1name="No QC", QC2name="QC on variant and sumstat", QC3name="QC on all", boot_num=50)

##UKBB
plot_R2_facut_allsumstat(kunkle_adsp,kunkle_adsp_qc,kunkle_UKBB_qc, bellenguez_adsp,bellenguez_adsp_qc, bellenguez_UKBB_qc, wightman_adsp, wightman_adsp_qc,wightman_UKBB_qc,
                         QC1name="No QC",  QC2name="QC on all", QC3name="QC, UKBB variants only",boot_num=50)


plot_auc_facut_all_sumstat(kunkle_adsp,kunkle_adsp_qc,kunkle_UKBB_qc, bellenguez_adsp,bellenguez_adsp_qc, bellenguez_UKBB_qc, wightman_adsp, wightman_adsp_qc,wightman_UKBB_qc,
                         QC1name="No QC",  QC2name="QC on all", QC3name="QC, UKBB variants only")


plot_auc_facut_all_sumstat(kunkle_adsp,kunkle_qc_variant_sumstat,kunkle_adsp_qc, bellenguez_adsp,  bellenguez_qc_on_variant_sumstat,bellenguez_adsp_qc, wightman_adsp,wightman_qc_variant_sumstat, wightman_adsp_qc,
                           QC1name="No QC", QC2name="QC on variant and sumstat", QC3name="QC on all")


## R2
## polyfun ------
mod1 <- glm(Diagnosis ~ Sex + Age + PRS, data= extract_afr(bellenguez_polypred_new), family=binomial)
RsqGLM(mod1, plot=FALSE)$Nagelkerke

#RsqGLM(glm(Diagnosis~PRS1+Sex+Age,family = 'binomial', data = bellenguez_polypred_new), plot=FALSE)$Nagelkerke
RsqGLM(glm(Diagnosis~PRS+Sex+Age,family = 'binomial', data = bellenguez_polypred_new), plot=FALSE)$Nagelkerke


plot_ethnic_R2(kunkle_polypred, "kunkle_polypred",  polypred_col, 50, "max_snp_per LD = ")
plot_ethnic_R2(wightman_polypred, "wightman_polypred",  polypred_col, 50, "max_snp_per LD = ")
plot_ethnic_R2(bellenguez_fixed_0224, "bellenguez_polypred",  polypred_col, 50, "max_snp_per LD = ")

## bellenguez----
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




