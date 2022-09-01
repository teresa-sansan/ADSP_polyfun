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
library(rsq)
#library(MESS)
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
wightman_susie_max10_polypred<- pre_process('/gpfs/commons/home/tlin/output/prs/wightman/susie_new_beta_max_snp_10_polypred.tsv')

###  polyfun-Pred ----
kunkle_polypred <- pre_process('/gpfs/commons/home/tlin/output/prs/polypred/kunkle/new_plink_polypred.tsv')
kunkle_susie <- pre_process('/gpfs/commons/home/tlin/output/prs/polypred/kunkle/new_plink_susie.tsv') ## 0811 add
kunkle_bl <- pre_process('/gpfs/commons/home/tlin/output/prs/polypred/kunkle/new_plink_bl_polypred.tsv')
kunkle_susie <- pre_process("/gpfs/commons/home/tlin/output/prs/polypred/kunkle/susie.prs.tsv")

bellenguez_polypred <- pre_process('/gpfs/commons/home/tlin/output/prs/polypred/bellenguez/old_plink_polypred.tsv')
wightman_polypred <- pre_process('/gpfs/commons/home/tlin/output/prs/polypred/wightman/fixed_0224.prs.tsv')

## new plink ------
bellenguez_adsp <- pre_process('/gpfs/commons/home/tlin/output/prs/new_plink/bellenguez/bellenguez_ADSP.tsv')
wightman_adsp <- pre_process('/gpfs/commons/home/tlin/output/prs/new_plink/wightman/wightman_ADSP.tsv')

bellenguez_adsp_qc <- pre_process('/gpfs/commons/home/tlin/output/prs/new_plink/bellenguez/bellenguez_ADSP_qc_all.tsv')
#wightman_adsp_qc <- pre_process('/gpfs/commons/home/tlin/output/prs/new_plink/wightman/wightman_ADSP_qc_all.tsv')
wightman_adsp_qc <- pre_process('/gpfs/commons/home/tlin/output/prs/new_plink/wightman/fixed_beta/wightman_ADSP_qc.tsv')

kunkle_UKBB_qc <- pre_process('/gpfs/commons/home/tlin/output/prs/new_plink/kunkle/kunkle_ADSP_UKBB_qc.tsv')
bellenguez_UKBB_qc <- pre_process('/gpfs/commons/home/tlin/output/prs/new_plink/bellenguez/bellenguez_ADSP_UKBB_qc.tsv')
wightman_UKBB_qc <- pre_process('/gpfs/commons/home/tlin/output/prs/new_plink/wightman/wightman_ADSP_UKBB_qc.tsv')

### susie -----
bellenguez_susie <- pre_process("/gpfs/commons/home/tlin/output/prs/polypred/bellenguez/susie_max_snp_10_polypred.tsv")


## polyfun_beta_pt -----
kunkle_polyfun_pT = pre_process('/gpfs/commons/home/tlin/output/prs/new_plink/kunkle/kunkle_polyfun_beta.tsv')
bellenguez_polyfun_pT = pre_process('/gpfs/commons/home/tlin/output/prs/new_plink/bellenguez/bellenguez_polyfun_beta.tsv')
wightman_polyfun_pT = pre_process('/gpfs/commons/home/tlin/output/prs/new_plink/wightman/wightman_polyfun_beta.tsv')

kunkle_polyfun_plink_no_cpT = pre_process('/gpfs/commons/home/tlin/output/prs/new_plink/kunkle/kunkle_polyfun_beta_noclump.tsv')
bellenguez_polyfun_plink_no_cpT = pre_process('/gpfs/commons/home/tlin/output/prs/new_plink/bellenguez/bellenguez_polyfun_beta_noclump.tsv')
wightman_polyfun_plink_no_cpT= pre_process('/gpfs/commons/home/tlin/output/prs/new_plink/wightman/wightman_polyfun_beta_noclump.tsv')

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
## Where method = how many causal SNPs are assumed in the LD block/ pT threshold.
plot_facet_dist <- function(df, title){
  ggplot(pivot_df(df), aes(x =PRS, group=Diagnosis_state,fill = Diagnosis_state, color = Diagnosis_state))+
    geom_density(alpha = 0.01) +
    facet_grid(method~final_population)+theme_bw()+ggtitle(title)
} 

## plot race/method separated dist plots
plot_facet_dist(kunkle_polypred, "kunkle_polypred")
plot_facet_dist(kunkle_new_beta, "kunkle_new_beta")
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

## Set column name & extract race ------
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
col_roc_polypred <- list("PRS1","PRS3","PRS5","PRS7","PRS10")
col_roc_polypred3 <- list("PRS1","PRS5","PRS10")

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
    #set.seed(10)
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
plot_ethnic_roc <- function(df, title, col,boot=TRUE,boot_num=5, plot=FALSE){
  if (boot != TRUE){ ## no bootstraping, use the original roc_result
    EUR = roc_result(extract_eur(df), title = paste(title, ", EUR") , column_for_roc = col, plot=plot)
    AFR = roc_result(extract_afr(df), title = paste(title, ", AFR") , column_for_roc = col, plot=plot)
    AMR = roc_result(extract_amr(df), title = paste(title, ", AMR") , column_for_roc = col, plot=plot)
    if (plot != FALSE){
      plot_grid(EUR, AFR, AMR,ncol = 1, nrow = 3)
      return() }
  } else {
    print(paste("boot time = ",boot_num))
    EUR = roc_result_boot(extract_eur(df), column_for_roc = col, boot_num = boot_num)
    AFR = roc_result_boot(extract_afr(df), column_for_roc = col, boot_num = boot_num)
    AMR = roc_result_boot(extract_amr(df), column_for_roc = col, boot_num = boot_num)
  }
  
  df <- data.frame(matrix(ncol = 0, nrow = length(col)*3))
  df$PRS =  rep(unlist(col),3)
  df$ethnicity = c(rep("EUR",length(col)),rep("AFR",length(col)),rep("AMR",length(col)))
  df$ethnicity <- factor(df$ethnicity,      # Reordering group factor levels
                         levels = c("EUR","AFR","AMR"))
  if (boot != TRUE){
    df$auc = append(append(unlist(EUR),unlist(AFR)),unlist(AMR))
  }
  else{
    df$auc = append(append(unlist(EUR$boot_mean),unlist(AFR$boot_mean)),unlist(AMR$boot_mean))
    df$boot_CI_lower = append(append(unlist(EUR$boot_CI_lower),unlist(AFR$boot_CI_lower)),unlist(AMR$boot_CI_lower))
    df$boot_CI_upper = append(append(unlist(EUR$boot_CI_upper),unlist(AFR$boot_CI_upper)),unlist(AMR$boot_CI_upper))
    df[c("auc","boot_CI_lower","boot_CI_upper")] <- sapply(df[c("auc","boot_CI_lower","boot_CI_upper")],as.numeric)
  }
  return(df)
}


plot_auc_facet_all_sumstat <- function(s1_qc1, s1_qc2, s1_qc3,s2_qc1, s2_qc2, s2_qc3,s3_qc1, s3_qc2, s3_qc3, col = col_roc_E5, QC1name="QC_on_base", QC2name="QC_on_base+variants",QC3name="no_qc",boot=TRUE, boot_num=5, legendname = FALSE){
  a = plot_ethnic_roc_facet(s1_qc1, s1_qc2, s1_qc3, col,'Kunkle',QC1name, QC2name,QC3name,boot=boot,boot_num=boot_num,legendname =legendname)
  b = plot_ethnic_roc_facet(s2_qc1, s2_qc2, s2_qc3, col,'Bellenguez',QC1name, QC2name,QC3name,boot=boot,boot_num=boot_num,legendname =legendname)
  c = plot_ethnic_roc_facet(s3_qc1, s3_qc2, s3_qc3, col,'Wightman',QC1name, QC2name,QC3name,boot=boot,boot_num=boot_num,legendname =legendname)
  
  print(a)
  prow <- plot_grid(a+ theme(legend.position="none"), 
                    b+ theme(legend.position="none",axis.text.y = element_blank())+ ylab(NULL), 
                    c+ theme(legend.position="none",axis.text.y = element_blank())+ylab(NULL),
                    ncol = 3, nrow = 1)
   legend <- get_legend(
     #c+ guides(legendname,color = guide_legend(nrow = 1)) +
     c +guides(col=guide_legend(legendname, nrow=1))+
       theme(legend.position = "bottom")
   )
   plot_grid(prow, legend, ncol=1,rel_heights=c(3,.4))
   
}

## fixed legend
## The boolean APOE operator is to see whether you want to plot the effect of only using 5(6, depending on QC or not) APOE allele to make prediction.  
plot_ethnic_roc_facet <- function(QC1, QC2, QC3, col=col_roc_E5, title=' ',QC1name="QC_on_base", QC2name="QC_on_base+variants",QC3name="no_qc",boot=TRUE, 
                                  boot_num=50, legendname=FALSE, APOE=FALSE){
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
  print(all)
  #all %>% 
  #  filter(PRS !="PRS-001" & PRS != "PRS_1") 
  if ('PRS_5' %in% col){
    all$PRS = str_replace(all$PRS, 'PRS_', 'p = 0.')
    all[all$PRS=='p = 0.e5',]$PRS =" p = 1e-5"
    
  }
  if ('PRS10' %in% col){
    all$PRS = str_replace(all$PRS, 'PRS', 'max snp ')
    
  }
  
  plot <- ggplot(data = all, aes(x=auc, y = PRS, color = qc_status))+
    geom_point(size=3,alpha=0.7,position = position_dodge(width = 0.7))+
    facet_wrap(~ethnicity, ncol=1)+
    xlab('AUC')+ ggtitle(title)+xlim(0.45, 0.75)+
    theme_bw()
  
  if (boot == TRUE){  
    plot <- plot + geom_errorbar(aes(xmin=boot_CI_lower, xmax=boot_CI_upper),position=position_dodge(width=0.7), width=.1,alpha=0.5,show.legend = FALSE) 
  }## plot error bar or not
  if(legendname != FALSE){ 
    plot = plot +guides(col=guide_legend(legendname))
  } ## change legend name or not (depends on usage)
  if(class(APOE) != "logical"){
    APOE_ci = plot_ethnic_roc(APOE, col="PRS", plot='F',boot=boot, boot_num=boot_num)
    #print(APOE_ci)
    plot = plot +  geom_vline(data = APOE_ci, aes(xintercept = auc,color='only APOE'),color="black",lwd=0.5, alpha=0.7,linetype="dashed")  ## move color into aes will generate legend automatically.
  } ## if we want to plot APOE (line) ## should do bootstrap automatically if it is specified 
  return(plot)
  #return(all)
}

## R2 functions ----
## calcualte pseudo rsquare 


rsq_formula <- function(formula, data, indices=FALSE) {    
  if (indices !=FALSE){
    d <- data[indices,] # allows boot to select sample
  }
  else 
    d = data 
  fit <- glm(formula = formula,family = 'binomial', data = d)
  #return(RsqGLM(fit, plot = FALSE)$Nagelkerke)
  return(fit)
} ## this is for boostrapping


log_reg_try_partial <- function(df,prs, boot_num = FALSE){
  mod_full <- glm(Diagnosis ~ Sex + Age + PRS5, family = 'binomial', data = bellenguez_fixed_0224) ## whole model
  mod_reduced <- glm(Diagnosis ~ PRS5, family = 'binomial', data = bellenguez_fixed_0224)  ## only PRS
  print(rsq.partial(mod_full, mod_reduced, type='n'))
  
}

log_reg_try_partial(bellenguez_fixed_0224, col_roc_polypred)




log_reg <- function(df,prs,main_title, plot = TRUE, legend="PRS_", boot_num = FALSE, replace="pT = 0."){
  set.seed(50)
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





print(la$partial.rsq)

## plot facet 
## plot ethnic facut plot for one sumstat
plot_ethnic_R2_facet <- function(QC1, QC2, QC3=FALSE, col=col_roc_E5, boot_num=50, title=' ',QC1name="QC_on_base", QC2name="QC_on_base+variants",QC3name="no_qc", legendname=FALSE){
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
  #all %>% 
  #  filter(PRS !="PRS-001" & PRS != "PRS_1") 
  if ('PRS_5' %in% col){
    all$PRS = str_replace(all$PRS, 'PRS_', 'p = 0.')
    all[all$PRS=='p = 0.e5',]$PRS =" p = 1e-5"
    
  }
  if ('PRS10' %in% col){
    all$PRS = str_replace(all$PRS, 'PRS', 'max snp ')
    
  }
  all$boot_mean = all$boot_mean * 100
  plot <- ggplot(data = all, aes(x= boot_mean, y = PRS, color = qc_status))+
    geom_point(size=3, alpha=0.75,position = position_dodge(width = 0.7))+
    geom_pointrange(aes(xmin=boot_CI_lower*100, xmax=boot_CI_upper*100), linetype="dotted",position=position_dodge(width=0.7),show.legend = FALSE) +
    facet_wrap(~ethnicity, ncol=1)+
    xlab('R squared (%)')+ ggtitle(title)+xlim(-1, 25)+
    theme_bw()
  if(legendname != FALSE){  ## testout new legend title
    plot = plot +guides(col=guide_legend(legendname))
  }
  return(plot)
}

## plot ethnic facut plot for multiple sumstat
plot_R2_facet_allsumstat<- function(s1_qc1, s1_qc2, s1_qc3,s2_qc1, s2_qc2, s2_qc3,s3_qc1, s3_qc2, s3_qc3,col = col_roc_E5, QC1name="QC_on_base", QC2name="QC_on_base+variants",QC3name="no_qc",boot_num=50, legendname=legendname){
  a = plot_ethnic_R2_facut(s1_qc1, s1_qc2, s1_qc3, col, boot_num, 'Kunkle', QC1name, QC2name,QC3name,legendname=legendname)
  b = plot_ethnic_R2_facut(s2_qc1, s2_qc2, s2_qc3, col, boot_num, 'Bellenguez',QC1name, QC2name,QC3name,legendname=legendname)
  c = plot_ethnic_R2_facut(s3_qc1, s3_qc2, s3_qc3, col, boot_num,'Wightman', QC1name, QC2name,QC3name,legendname=legendname)
  prow <- plot_grid(a+ theme(legend.position="none"), 
                    b+ theme(legend.position="none",axis.text.y = element_blank())+ ylab(NULL), 
                    c+ theme(legend.position="none",axis.text.y = element_blank())+ylab(NULL),
                    ncol = 3, nrow = 1)
  legend_b <- get_legend(
    c+ guides(color = guide_legend(legendname, nrow = 1)) +
      theme(legend.position = "bottom")
  )
  plot_grid(prow, legend_b, ncol=1,rel_heights=c(3,.4))
  
}

## try out partial R2 ----



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
plot_ethnic_roc_facet(bellenguez_cT0224,bellenguez_adsp,bellenguez_adsp_qc,col_roc_E5,title='bellenguez',
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

plot_ethnic_roc_facet(kunkle_qc_all, kunkle_withoutAPOE_qc,title='kunkle (QC)',
                      QC1name = "kunkle", QC2name = "kunkle_no_APOE", 
                      legendname = "APOE status", boot_num=2)

kunkle_pt_qc_plot_maf=roc_result(kunkle_qc_maf,title="QC, MAF >0 0.1%", column_for_roc = col_roc_E5)
##kunkle_pt_qc_plot_base=roc_result(kunkle_qc_base, title="kunkle_qc_summary_stats",column_for_roc = col_roc_E5)
plot_grid(kunkle_pt_qc_plot,kunkle_pt_qc_plot_maf,ncol = 2, nrow = 2)


## APOE
plot_ethnic_roc(kunkle_withoutAPOE, "kunkle remove APOE region", col_roc_E5)
col_roc <- list("only2SNP","no2SNP","no2SNP_qc")

kunkle_APOE_roc <- roc(Diagnosis ~ only2SNP+no2SNP+no2SNP_qc, data = kunkle_APOE)
auc(roc(Diagnosis~only2SNP, data = kunkle_APOE)) ##0.621
auc(roc(Diagnosis~only2SNP, data = extract_eur(kunkle_APOE))) ##0.6187

for (i in list(extract_eur(kunkle_APOE),extract_afr(kunkle_APOE),extract_amr(kunkle_APOE))){
  print(auc(roc(Diagnosis~no2SNP, data = i, quiet=T))[1])
}

plot_ethnic_R2_facet(kunkle_adsp_qc, kunkle_withoutAPOE_qc, FALSE, boot_num = 50,title='Kunkle (QC)',
                     QC1name='withAPOE', QC2name='without APOE', legendname='APOE status')


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


## polyfun_AUC----
##kunkle
plot_ethnic_roc_facet(kunkle_bl, kunkle_polypred,FALSE,col =col_roc_polypred3, title= 'kunkle (new_plink)',
                      QC1name = 'Baseline annotations', QC2name = 'All annotations',
                      legendname = 'Annotation'
)


plot_ethnic_roc_facet(kunkle_susie, kunkle_bl, kunkle_polypred, col = col_roc_polypred3, title= 'kunkle',
                      QC1name = 'SuSiE', QC2name = 'PolyFun (BL)',QC3name='PolyFun (All anno)',
                      legendname = 'Annotation'
)


plot_ethnic_R2_facet(kunkle_susie, kunkle_bl, kunkle_polypred, col =col_roc_polypred3, title= 'kunkle (new_plink)',
                     QC1name = 'no annotations(SuSiE)', QC2name = 'Baseline annotations', QC3name = 'All annotations',
                     legendname = 'Annotation'
)

plot_ethnic_roc_facet(bellenguez_susie, bellenguez_polypred,FALSE,col =col_roc_polypred3, title= 'bellenguez (old_plink)',
                      QC1name = 'SuSiE', QC2name = 'Polyfun',
                      legendname = 'Annotation'
)


plot_ethnic_roc_facet(bellenguez_susie, bellenguez_polypred,FALSE,col =col_roc_polypred3, title= 'bellenguez (old_plink)',
                      QC1name = 'SuSiE', QC2name = 'Polyfun',
                      legendname = 'Annotation'
)
wightman_polypred['PRS']= wightman_polypred['PRS10']

plot_ethnic_R2_facet(wightman_susie_max10, wightman_polypred,FALSE, col ='PRS', title= 'wightman',
                     QC1name = 'SuSiE', QC2name = 'Polyfun',
                     legendname = 'Annotation'
)

plot_ethnic_roc_facet(bellenguez_polypred_new, bellenguez_susie, FALSE, boot_num = 50, col='PRS',
                      title='bellenguez', QC1name = 'PolyFun', QC2name = 'SuSiE', legendname = 'Tools')

plot_ethnic_R2_facet(bellenguez_polypred_new, bellenguez_susie, FALSE, boot_num = 50, col='PRS',
                      title='bellenguez', QC1name = 'PolyFun', QC2name = 'SuSiE', legendname = 'Tools')

## polyfun beta vs sumstat beta ---
plot_auc_facet_all_sumstat(kunkle_adsp_qc, kunkle_polyfun_pT, kunkle_polyfun_plink_no_cpT, 
                           bellenguez_adsp_qc, bellenguez_polyfun_pT,bellenguez_polyfun_plink_no_cpT,
                           wightman_adsp_qc, wightman_polyfun_pT, wightman_polyfun_plink_no_cpT,
                           QC1name = 'Summary Stat', QC2name = 'PolyFun (c+pT)', QC3name = 'PolyFun', legendname = 'Effect size', boot_num=50
                           )

plot_R2_facet_allsumstat(kunkle_adsp_qc, kunkle_polyfun_pT, kunkle_polyfun_plink_no_cpT, 
                           bellenguez_adsp_qc, bellenguez_polyfun_pT,bellenguez_polyfun_plink_no_cpT,
                           wightman_adsp_qc, wightman_polyfun_pT, wightman_polyfun_plink_no_cpT,
                           QC1name = 'Summary Stat', QC2name = 'PolyFun (c+pT)', QC3name = 'PolyFun',legendname = 'Effect size', boot_num=50
)

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

plot_ethnic_roc_facet(kunkle_cT,kunkle_withoutAPOE,kunkle_withoutAPOE_qc,col_roc_E5, title='kunkle',
                      QC1name= "with APOE", QC2name= "without APOE", QC3name="without APOE, QC")

plot_auc_facet_all_sumstat(kunkle_qc_base, kunkle_qc_variant_sumstat, kunkle_cT,bellenguez_qc_on_base, bellenguez_qc_on_variant_sumstat,bellenguez_cT0224, 
                           wightman_qc_base, wightman_qc_variant_sumstat, wightman_cT) 

#plot_auc_facut_all_sumstat(kunkle_adsp,kunkle_qc_variant_sumstat,kunkle_adsp_qc, bellenguez_adsp,  bellenguez_qc_on_variant_sumstat,bellenguez_adsp_qc, wightman_adsp,wightman_qc_variant_sumstat, wightman_adsp_qc)
plot_R2_facet_allsumstat(kunkle_adsp,kunkle_qc_variant_sumstat,kunkle_adsp_qc, bellenguez_adsp,  bellenguez_qc_on_variant_sumstat,bellenguez_adsp_qc, wightman_adsp,wightman_qc_variant_sumstat, wightman_adsp_qc,
                         QC1name="No QC", QC2name="QC on variant and sumstat", QC3name="QC on all", boot_num=50)

##UKBB
plot_R2_facet_allsumstat(kunkle_adsp,kunkle_adsp_qc,kunkle_UKBB_qc, bellenguez_adsp,bellenguez_adsp_qc, bellenguez_UKBB_qc, wightman_adsp, wightman_adsp_qc,wightman_UKBB_qc,
                         QC1name="No QC",  QC2name="QC on all", QC3name="QC, UKBB variants only",boot_num=50)


plot_auc_facet_all_sumstat(kunkle_adsp,kunkle_adsp_qc,kunkle_UKBB_qc, bellenguez_adsp,bellenguez_adsp_qc, bellenguez_UKBB_qc, wightman_adsp, wightman_adsp_qc,wightman_UKBB_qc,
                           QC1name="No QC",  QC2name="QC on all", QC3name="QC, UKBB variants only")


plot_auc_facet_all_sumstat(kunkle_adsp,kunkle_qc_variant_sumstat,kunkle_adsp_qc, bellenguez_adsp,  bellenguez_qc_on_variant_sumstat,bellenguez_adsp_qc, wightman_adsp,wightman_qc_variant_sumstat, wightman_adsp_qc,
                           QC1name="No QC", QC2name="QC on variant and sumstat", QC3name="QC on all")

## R2
## polyfun ------

## bellenguez----
#df,prs,main_title, plot = TRUE, legend="PRS_", boot_num = FALSE, replace="pT = 0."
plot_R2_cT(bellenguez_cT, bellenguez_qc, bellenguez_qc_on_base, bellenguez_qc_on_target, 'bellenguez_updateRSID', col_roc, plot=TRUE, boot_num = 50)
plot_R2_cT(bellenguez_cT0224, bellenguez_qc0224, bellenguez_qc_on_base0224, bellenguez_qc_on_target0224, 'bellenguez_fixed0224', col_roc_E5, boot_num=50 )

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

## wightman ------
plot_ethnic_R2(wightman_cT, col_roc_E5, 'Wightman', 50)


##susie ------
plot_ethnic_R2(susie, col_roc_polypred, "bellenguez_susie",50)
plot_ethnic_R2(kunkle_susie, col_roc_polypred, "kunkle_susie",50)

RsqGLM(glm(Diagnosis~PRS10+Sex+Age,family = 'binomial', data = kunkle_susie), plot=FALSE)$Nagelkerke ## 0.14
RsqGLM(glm(Diagnosis~PRS1+Sex+Age,family = 'binomial', data = kunkle_susie), plot=FALSE)$Nagelkerke ##0.139

##others -----
sbayesr_log <-  glm(Diagnosis~PRS+Sex+Age,family = 'binomial', data = sbayesR)
prsice_log <- glm(Diagnosis~PRS+Sex+Age,family = 'binomial', data = PRSice)

# PRS_dist density
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

plot_density <- function(df1, df2, df3, df4, sumstat_name, Pthreshold, df1_name = '', df2_name = '', df3_name = '', df4_name = '', eur=FALSE){
  #par(mfrow=c(2,2),xpd="NA")
  if(eur==TRUE){
    print("plot EUR subset")
    df = extract_eur(df1)
    df_qc = extract_eur(df2)
    df_qc_base = extract_eur(df3)
    df_qc_target = extract_eur(df4)
  }
  PRS_density(df1,  paste(sumstat_name, ", ",  df1_name), Pthreshold)
  PRS_density(df2,  paste(sumstat_name, ", ", df2_name), Pthreshold)
  PRS_density(df3,  paste(sumstat_name, ", ", df3_name),Pthreshold)
  PRS_density(df4, paste(sumstat_name, ", ", df4_name),Pthreshold)
}

plot_density(kunkle_cT, kunkle_qc, kunkle_qc_base, kunkle_qc_target, "kunkle_eur", "PRS_05", eur=TRUE)
legend(-0.32,69, inset=.02, title="AD diagnosis",
       c("Case","Control"), fill=c("red","black"), horiz=TRUE, cex=0.8,xpd="NA")
plot_density(wightman_qc,wightman_qc_target,wightman_qc_variant,wightman_qc_individual, 'wightman','PRS_5',
             df1_name = 'qc',df2_name = 'qc target',df3_name = 'qc variant',df4_name = 'qc individual')

legend("topleft", inset = c(- 0.4),                  
       legend = c("Case","Control"), col = 1:2,fill=c("red","black"), horiz=TRUE, cex=0.8, title= 'AD diagnosis')

### plot density/age/case-----
DATA = kunkle_withoutAPOE_qc %>%
  subset(final_population != "SAS",) %>%
  subset(final_population != "EAS",) %>%
  subset(final_population != "UNKNOWN",)


DATA = kunkle_adsp_qc %>%
  subset(final_population != "SAS",) %>%
  subset(final_population != "EAS",) %>%
  subset(final_population != "UNKNOWN",)


DATA$Diagnosis = as.factor(DATA$Diagnosis)
DATA$Diagnosis_state=ifelse(DATA$Diagnosis == 0, "control","case") 

## function provide by https://stackoverflow.com/questions/35717353/split-violin-plot-with-ggplot2 
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
                             data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                             grp <- data[1, "group"]
                             newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                             newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                             newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
                             
                             if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                               stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                                         1))
                               quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                               aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
                               aesthetics$alpha <- rep(1, nrow(quantiles))
                               both <- cbind(quantiles, aesthetics)
                               quantile_grob <- GeomPath$draw_panel(both, ...)
                               ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                             }
                             else {
                               ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                             }
                           })


geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
        position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}

ggplot(data=DATA, aes(x=factor(final_population, level = c('EUR','AFR','AMR')), y=Age, fill=Diagnosis_state))+
  geom_split_violin(alpha = 0.3)+geom_boxplot(width=0.4, alpha=0.8,show.legend = FALSE)+ xlab('Ethnicity')+  scale_fill_discrete(name ="AD diagnosis")+coord_flip()+
  theme_bw()


## draw APOE analysis (target file) -----

roc(kunkle_APOE[["Diagnosis"]],kunkle_APOE$PRS, quiet=T)

## AUC = as.numeric(auc(Diagnosis~PRS))

APOE_subtype = function(df){
  new_df = df %>% 
    subset(final_population != 'SAS' & final_population != 'EAS' & final_population != 'UNKNOWN') %>% 
    group_by(APOE,Diagnosis, final_population) %>%
    summarise(PRS_mean = mean(PRS), PRS_std =sd(PRS), Age_mean = mean(Age))
  new_df$Diagnosis = as.factor(new_df$Diagnosis)
  levels(new_df$Diagnosis) <- c("Control", "Case") 
  return(new_df)}

APOE = APOE_subtype(kunkle_withoutAPOE_qc)
APOE = APOE_subtype(kunkle_adsp_qc)

ggplot(data = APOE, aes(y  = PRS_mean, x = Age_mean, color = APOE, shape = Diagnosis)) + 
  geom_point(size=3, alpha = 0.8) + facet_wrap(~factor(final_population, levels = c("EUR","AFR","AMR")), ncol=1)+
  theme_bw()+ggtitle("kunkle without APOE (QC)")

ggplot(data = APOE, aes(y  = PRS_mean, x = APOE, fill=Diagnosis)) + 
  geom_bar(stat="identity", position=position_dodge(), alpha=0.8) + facet_wrap(~factor(final_population, levels = c("EUR","AFR","AMR")), ncol=1)+
  scale_fill_brewer(palette="Blues")+
  theme_bw()


