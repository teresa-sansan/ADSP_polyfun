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
library(devtools)
library(terra)
library(modEvA)  # psuedo rsquare
library(boot)
library(gridGraphics)
library(UpSetR)

pre_process <- function(df, file=FALSE){
  if(file==FALSE){
    df<- read.csv(df,sep = '\t', header=T,fill = T)
  }
  if("age_covariate" %in% colnames(df))
  {
    print("change column name first...")
    colnames(df)[which(names(df) == 'age_covariate')] <- 'Age'
    colnames(df)[which(names(df) == 'AD_status_final')] <- 'Diagnosis'
  }
  df$Age <- as.numeric(as.character(df$Age))
  print(paste("original=",  dim(df)[1], "rows"))
  #df <- df %>%
  #  filter(Diagnosis != -1 & Age >= 65)
  df <- df[df$Age >= 65,]
  print(paste("filtered=",  dim(df)[1], "rows"))
  df$Diagnosis <- as.numeric(df$Diagnosis)
  print(colnames(df))
  if("PRS_0.1" %in% colnames(df))
  {
    print("change PRS column name")
    colnames(df) = gsub("_0\\.", "_", colnames(df))
  }
  return(df)
} ## remove the sample younger than 65 || have no diagnosis || rename col

## load PRS data ----
### pT -----

###  36k ----
wightman_new_anno_0824_only_ml_ibd<-  pre_process('/gpfs/commons/home/tlin/output/prs/polypred/new_anno_0824/wightman.no_ml.ibd36k.tsv')
wightman_new_anno_0824_all_ibd <-  pre_process('/gpfs/commons/home/tlin/output/prs/polypred/new_anno_0824/wightman.all.ibd36k.tsv')
wightman_new_anno_0824_no_ml_ibd <-  pre_process('/gpfs/commons/home/tlin/output/prs/polypred/new_anno_0824/wightman.no_ml.ibd36k.tsv')
wightman_new_anno_0824_bl_ibd <-  pre_process('/gpfs/commons/home/tlin/output/prs/polypred/new_anno_0824/wightman.bl.ibd36k.tsv')
wightman_new_anno_0824_susie_ibd <- pre_process('/gpfs/commons/home/tlin/output/prs/polypred/new_anno_0824/wightman.susie.adj_beta_ibd36k.tsv')
wightman_allanno_ibd <- pre_process('/gpfs/commons/home/tlin/output/wightman/new_anno_0203/wightman.update_all+enformer.ibd36k.tsv')

bellenguez_new_anno_0824_all_ibd <-  pre_process('/gpfs/commons/home/tlin/output/prs/polypred/new_anno_0824/bellenguez_all_ibd36k.tsv')
bellenguez_new_anno_0824_no_ml_ibd <-  pre_process('/gpfs/commons/home/tlin/output/prs/polypred/new_anno_0824/bellenguez_no_ml_ibd36k.tsv')
bellenguez_new_anno_0824_only_ml_ibd <-  pre_process('/gpfs/commons/home/tlin/output/prs/polypred/new_anno_0824/bellenguez_only_ml_ibd36k.tsv')
bellenguez_new_anno_0824_bl_ibd <-  pre_process('/gpfs/commons/home/tlin/output/prs/polypred/new_anno_0824/bellenguez_bl_ibd36k.tsv')
bellenguez_new_anno_0824_susie <-  pre_process('/gpfs/commons/home/tlin/output/prs/polypred/new_anno_0824/bellenguez_susie_ibd36k.tsv')

## new_anno 0318_24
kunkle_0318 <- pre_process('/gpfs/commons/home/tlin/output/prs/new_anno_0318_24/kunkle/kunkle.tsv')
wightman_0318 <- pre_process('/gpfs/commons/home/tlin/output/prs/new_anno_0318_24/wightman/wightman.tsv')
bellenguez_0318 <- pre_process('/gpfs/commons/home/tlin/output/prs/new_anno_0318_24/bellenguez/bellenguez.tsv')
bellenguez_adsp_0318 <- pre_process('/gpfs/commons/home/tlin/output/prs/new_anno_0318_24/bellenguez_adsp_reference/bellenguez_adsp_reference.tsv')


## new_anno_0318_24 with thres
bellenguez_adsp_0318_no_thres <- pre_process('/gpfs/commons/home/tlin/output/prs/new_anno_0318_24/bellenguez_adsp_reference/thres/prs_no_pip_thres.tsv')
bellenguez_adsp_0318_pip01 <- pre_process('/gpfs/commons/home/tlin/output/prs/new_anno_0318_24/bellenguez_adsp_reference/thres/prs_pip0.1.tsv')
bellenguez_adsp_0318_pip025 <- pre_process('/gpfs/commons/home/tlin/output/prs/new_anno_0318_24/bellenguez_adsp_reference/thres/prs_pip0.25.tsv')


### Other PRS method -----
PRSice <- pre_process("/gpfs/commons/home/tlin/output/prs/PRSice_pheno.tsv")
sbayesR = pre_process("/gpfs/commons/home/tlin/output/prs/sbayesR.tsv")

### PRSCS ----
PRSCS <- pre_process('/gpfs/commons/home/tlin/output/wightman/prscs/original/prscs_17k.tsv')
polypred_PRSCS = pre_process('/gpfs/commons/home/tlin/output/wightman/prscs/all_anno/prscs_17k.tsv')
polypred_plink = pre_process('/gpfs/commons/home/tlin/output/wightman/prscs/all_anno/prs_plink_17k.tsv')
polypred_PRSCS_NOT0 <-  pre_process('/gpfs/commons/home/tlin/output/wightman/prscs/all_anno/PIP_not0/prscs_pipNOT0.tsv')
polypred_prscs_noEnformer <- pre_process('/gpfs/commons/home/tlin/output/wightman/prscs/all_except_enformer/prscs_pipNOT0.tsv')

##didnt use the BETA from polyfun (polyfun -> subset -> PRSCS -> plink)
polyfun_PRSCS_NOT0 <- pre_process('/gpfs/commons/home/tlin/output/wightman/prscs/all_anno/beta_sumstat/prscs_pipNOT0.tsv')

## other subset (subset[polyfun, PRSCS] -> plink )
PRSCS_subset <- pre_process("/gpfs/commons/home/tlin/output/wightman/prscs/original/subset_polyfun/prscs_pipNOT0.tsv")

#36k
PRSCS_36K <- pre_process('/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_preview/PRS_hg38/prscs/prs_36k.tsv')
PRSCS_36K <- merge(PRSCS_36K, lookup_table, by = "Race", all.x = TRUE)
PRSCS_36k_fixed_hispanic <- read.csv('/gpfs/commons/home/tlin/output/36k/wightman/PRSCS_fixed_hispanic.tsv', sep = '\t', header=T,fill = T)
PRSCS_36k_new_ancestry <- pre_process('/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_preview/PRS_hg38/prscs/prscs_new_ancestry.tsv')
PRSCS_36K_bellenguez<- pre_process('/gpfs/commons/home/tlin/output/prs/PRSCS/36k/bellenguez/prscs_36k.tsv')
PRSCS_36k_ibd_wightman <- pre_process('/gpfs/commons/home/tlin/output/prs/PRSCS/36k_ibd_adsp_fixed/wightman/prscs_wightman_ADSP_ibd_36k.tsv')

#Polypred (PRSCS_POLYFUN)
polypred <- pre_process('/gpfs/commons/home/tlin/output/wightman/new_anno_0203/update_all+enformer/finemap/polypred/w_prscs/polypred.prs.wpheno.tsv')


##CasioPR36K

inf_inter2_0219 <- pre_process('/gpfs/commons/home/tlin/pic/casioPR/simulation/prs/0219_inf_intersect2_whole_genome.tsv')
inf_0214<- pre_process('/gpfs/commons/home/tlin/pic/casioPR/simulation/prs/0214_inf_whole_genome.tsv')

##Extract diff age of controls -----
segregate_control_by_age <- function(df, age, balanced = TRUE){
  sprintf('Only keep controls whose age is > %d ... \n',age)
  set.seed(10)
  case = df[df$Diagnosis == 1,]
  case_num = dim(case)[1]
  control = df %>%
    filter(Diagnosis == 0) %>%
    filter(Age >= age)
  control_num = dim(control)[1]
  sprintf("Have %d cases and %d controls. \n",case_num,control_num)
  if(balanced == TRUE){
    sprintf("Balancing number of case and control ...\n")
    if(case_num > control_num){
      sprintf("Removing %d of case. ",case_num-control_num)
      sample_match = case[sample(nrow(case), size = control_num),]
      df_match = rbind(sample_match, control)
    }
    else{
      sprintf('Removing %d of control',control_num-case_num)
      sample_match = control[sample(nrow(control), size = case_num),]
      df_match = rbind(case, sample_match)
    }
    sprintf("Ending up with %d sample in total.", dim(df_match)[1])
    return(df_match)
  }
  print('returning not balanced result')
  return(rbind(case, control))
}

## Set prs col names ------
col_roc <- list("PRS_001","PRS_005","PRS_01","PRS_05","PRS_1","PRS_5")
col_roc_E5 <- list("PRS_e5","PRS_001","PRS_005","PRS_01","PRS_05","PRS_1","PRS_5")
col_roc_polypred <- list("PRS1","PRS3","PRS5","PRS7","PRS10")
col_roc_polypred3 <- list("PRS1","PRS5","PRS10")
col_roc_polypred_125<- list("PRS1","PRS2","PRS5")
col_roc_test_anno <-list ('PRS_bl','PRS_bl_omics','PRS_bl_omics_dl')
col_roc_anno <-list ('PRS_susie','PRS_bl','PRS_bl_omics','PRS_bl_omics_dl')

plot_ethnic_roc(PRSCS_36k_ibd_bellengeuz, col = list("PRS_1","PRS_5","PRS_01"))
plot_ethnic_roc(PRSCS_36k_ibd_wightman, col = list("PRS_1","PRS_5","PRS_01",'PRS_05'))

##result ----
# c+T (before and after qc)-----
## all
plot_auc_facet_all_sumstat(kunkle_adsp_no_apoe, kunkle_adsp_no_apoe_qc, FALSE, 
                           bellenguez_adsp, bellenguez_adsp_qc, FALSE,
                           wightman_UKBB, wightman_UKBB_qc, FALSE, 
                           jansen_adsp, jansen_adsp_qc, FALSE, legendname="QC status",
                           QC1name = "no qc",
                           QC2name = "qc", 
                           boot_num = 50)

## all without bellenguez
plot_auc_facet_all_sumstat(kunkle_adsp_no_apoe, kunkle_adsp_no_apoe_qc, FALSE, 
                           FALSE, FALSE, FALSE,
                           wightman_UKBB, wightman_UKBB_qc, FALSE, 
                           jansen_adsp, jansen_adsp_qc, FALSE,legendname='QC status',
                           QC1name = "no qc",
                           QC2name = "qc", 
                           boot_num = 50)

plot_ethnic_roc_facet(PRSCS, polypred_PRSCS, polypred_plink,title='wightman',
                      QC1name = "PRSCS", QC2name = "polyfun_prscs", 
                      legendname = "polypred_plink", boot_num=50)


## all partial R2
## no qc
for ( i in list(kunkle_adsp,kunkle_adsp_no_apoe, bellenguez_adsp,wightman_UKBB, jansen_adsp)){
  print("EUR")
  print(log_reg_partial(extract_eur(i), col_roc_E5))
  print("AFR")
  print(log_reg_partial(extract_afr(i), col_roc_E5))
  print("AMR")
  print(log_reg_partial(extract_amr(i), col_roc_E5))
  print("")
}

## qc
for ( i in list(kunkle_adsp_qc,kunkle_adsp_no_apoe_qc, bellenguez_adsp_qc,wightman_UKBB_qc, jansen_adsp_qc)){
  print("EUR")
  print(log_reg_partial(extract_eur(i), col_roc_E5))
  print("AFR")
  print(log_reg_partial(extract_afr(i), col_roc_E5))
  print("AMR")
  print(log_reg_partial(extract_amr(i), col_roc_E5))
  print("")
}

for ( i in list(wightman_bl, wightman_polypred, wightman_no_ml, wightman_all_except_enformer, wightman_glasslab, wightman_update_all)){
  print("EUR")
  print(log_reg_partial(extract_eur(i), col_roc_E5))
  print("AFR")
  print(log_reg_partial(extract_afr(i), col_roc_E5))
  print("AMR")
  print(log_reg_partial(extract_amr(i), col_roc_E5))
  print("")
}

## polyfun ------



#36k
plot_ethnic_roc_facet(bellenguez_new_anno_0824_susie, bellenguez_new_anno_0824_no_ml_ibd ,QC3=bellenguez_new_anno_0824_only_ml_ibd,QC4=bellenguez_new_anno_0824_all_ibd,
                      QC1name = "none", QC2name = "Omics",QC3name = 'DL', QC4name='Omics+DL',
                      col = list("PRS"),boot_num = 50, title='bellenguez 36k',legendname = 'annotation')




roc <- plot_ethnic_roc_facet(bellenguez_new_anno_0824_susie, bellenguez_new_anno_0824_no_ml_ibd ,QC3=bellenguez_new_anno_0824_only_ml_ibd,QC4=bellenguez_new_anno_0824_all_ibd,
                      QC1name = "none", QC2name = "Omics",QC3name = 'DL', QC4name='Omics+DL',
                      col = list("PRS"),boot_num = 50, title='bellenguez ',legendname = 'annotation')


plot <- ggplot(data = roc, aes(x=auc, y = qc_status, color = qc_status))+
  geom_point(size=3,alpha=0.7,position = position_dodge(width = 0.7))+
  facet_wrap(~ethnicity, ncol=1)+
  xlab('AUC')+ ggtitle('bellenguez')+xlim(0.47, 0.60)+
  theme_bw()


plot <- plot + geom_errorbar(aes(xmin=boot_CI_lower, xmax=boot_CI_upper),position=position_dodge(width=0.7), width=.1,alpha=0.5,show.legend = FALSE) 

plot_ethnic_roc_facet(wightman_new_anno_0824_susie_ibd, wightman_new_anno_0824_no_ml_ibd ,QC3=wightman_new_anno_0824_only_ml_ibd,QC4=wightman_new_anno_0824_all_ibd,
                      QC1name = "none", QC2name = "Omics",QC3name = 'DL', QC4name='Omics+DL',
                      col = list("PRS5"),boot_num = 50, title='wightman 36k, +age+sex',legendname = 'annotation')

plot_ethnic_roc_facet(wightman_new_anno_0824_susie_ibd, wightman_new_anno_0824_no_ml_ibd ,QC3=wightman_new_anno_0824_only_ml_ibd,QC4=wightman_new_anno_0824_all_ibd,
                      QC1name = "none", QC2name = "Omics",QC3name = 'DL', QC4name='Omics+DL',
                      col = list("PRS5"),boot_num = 50, title='wightman 36k',legendname = 'annotation')

## PRSCS ------
# Check the levels of the 'Diagnosis' variable
# Convert 'Diagnosis' from factor to numeric
polypred_PRSCS_NOT0$Diagnosis <- as.numeric(polypred_PRSCS_NOT0$Diagnosis)
polypred_prscs_noEnformer$Diagnosis <- as.numeric(polypred_prscs_noEnformer$Diagnosis)


## all, polyfun vs susie-----
plot_auc_facet_all_sumstat(kunkle_susie, kunkle_polypred, FALSE, 
                           FALSE, FALSE,FALSE,
                           wightman_susie, wightman_polypred, FALSE, 
                           jansen_susie, jansen_polypred, FALSE, legendname = "annotation",
                           col = col_roc_polypred3,
                           QC1name = "none",
                           QC2name = "all",
                           boot_num = 50)


## PRSCS vs Polyfun vs Polypred---
polypred_roc=plot_ethnic_roc(polypred, col='PRS', plot=FALSE, boot_num=50)
polypred_roc$PRS='PolyPred'
PRSCS_roc=plot_ethnic_roc(PRSCS, col='PRS_05', plot=FALSE, boot_num=50)
PRSCS_roc$PRS = "PRSCS"
polyfun_roc=plot_ethnic_roc(wightman_all_anno, col='PRS5', plot=FALSE, boot_num=50)
polyfun_roc$PRS = 'PolyFun'

df = rbind(polyfun_roc, polypred_roc, PRSCS_roc)
ggplot(data = df, aes(x=auc, y = PRS))+
  geom_point(size=2,alpha=0.9,position = position_dodge(width = 0.7), color='darkblue')+ylab('method')+
  xlab('AUC')+ ggtitle('Wightman')+xlim(0.3, 0.75)+theme_bw() +facet_wrap(~ethnicity, ncol=1)+
  geom_errorbar(aes(xmin=boot_CI_lower, xmax=boot_CI_upper),position=position_dodge(width=0.7), width=.1,alpha=0.5,color='darkblue',show.legend = FALSE) 


## partial r2
data_list <- list(kunkle_polypred, bellenguez_polypred, wightman_polypred, jansen_polypred)
lapply(data_list, apply_log_reg_partial)



## age segrgation

plot_control_age_roc_multi( plot_control_age_roc(bellenguez_new_anno_0824_susie, col = list("PRS"), plot=FALSE),
                            plot_control_age_roc(bellenguez_new_anno_0824_no_ml_ibd, col = list("PRS"), plot=FALSE),
                            plot_control_age_roc(bellenguez_new_anno_0824_all_ibd, col = list("PRS"), plot=FALSE),
                            names=list("none", "omics",'omics+DL'), title='Bellenguez 36k')




plot_partialR <- function(df, col = list('PRS')){
  eur = log_reg_partial_cov(extract_race(df,'EUR'), col )
  afr = log_reg_partial_cov(extract_race(df,'AFR'), col )
  amr = log_reg_partial_cov(extract_race(df,'AMR'), col )
  df_out = rbind(eur, afr,amr)
  colnames(df_out)[1] <- "ethnicity"
  df_out$ethnicity = c('EUR','AFR','AMR')
  return(df_out)
}

df1 = plot_partialR(bellenguez_new_anno_0824_susie)
df2 = plot_partialR(bellenguez_new_anno_0824_no_ml_ibd)
df3 = plot_partialR(bellenguez_new_anno_0824_only_ml_ibd)
df4 = plot_partialR(bellenguez_new_anno_0824_all_ibd)
df_out = rbind(df1,df2, df3,df4)
df_out['annotation'] = c('none','Omics','DL','Omics+DL')

df_out$annotation = factor(df_out$annotation, levels = unique( df_out$annotation ))
df_out$ethnicity = factor(df_out$ethnicity, levels = c('EUR','AFR','AMR'))
df_out_bellenguez= df_out
ggplot(data =df_out, aes(x = partial_R2_noage, y = annotation, color = annotation)) + 
  geom_point(size=2, alpha=0.9,position = position_dodge(width = 0.7)) + facet_wrap(~ethnicity, scales = "free_y", ncol=1)+
  labs( x = "partial R2(%)")+ggtitle('bellenguez36k')+
  theme_bw()

##wightman
df1 = plot_partialR(wightman_new_anno_0824_susie, col =list('PRS5'))
df2 = plot_partialR(wightman_new_anno_0824_no_ml_ibd, col =list('PRS5'))
df3 = plot_partialR(wightman_new_anno_0824_only_ml_ibd, col =list('PRS5'))
df4 = plot_partialR(wightman_new_anno_0824_all_ibd, col =list('PRS5'))
df_out = rbind(df1,df2, df3,df4)
df_out['annotation'] = c('none','Omics','DL','Omics+DL')

df_out$annotation = factor(df_out$annotation, levels = unique( df_out$annotation ))
df_out$ethnicity = factor(df_out$ethnicity, levels = c('EUR','AFR','AMR'))
df_out_wightman= df_out
ggplot(data =df_out_wightman, aes(x = partial_R2_noage, y = annotation, color = annotation)) + 
  geom_point(size=2, alpha=0.9,position = position_dodge(width = 0.7)) + facet_wrap(~ethnicity, scales = "free_y", ncol=1)+
  labs( x = "partial R2(%)")+ggtitle('wightman36k')+
  theme_bw()



## 
#36K_test_anno_chirag
plot_ethnic_roc_facet(bellenguez_0318 ,bellenguez_adsp_0318, QC3 =data.frame(),QC1name = "bellenguez_UKBB", QC2name = "bellenguez_adsp",
                      col = col_roc_test_anno,boot_num = 50, title='36k IBD test',legendname = 'sumstat')


plot_ethnic_roc(kunkle_0318, col=col_roc_test_anno, plot=TRUE,title='Kunkle') 
plot_ethnic_roc(wightman_0318, col=col_roc_test_anno, plot=TRUE,title='wightman') 
plot_ethnic_roc(bellenguez_0318, col=col_roc_test_anno, plot=TRUE,title='bellenguez') 

##adsp reference panel, with PIP thres
plot_ethnic_roc_facet(bellenguez_adsp_0318_no_thres ,bellenguez_adsp_0318_pip01, QC3 =bellenguez_adsp_0318_pip025,QC1name = "no thres", QC2name = "pip > 0.1",QC3name='pip > 0.25',
                      col = col_roc_anno,boot_num = 50, title='bellenguez adsp panel',legendname = 'pip thres')

plot_ethnic_roc(bellenguez_adsp_0318_no_thres, col=col_roc_anno, title='no thres', plot=TRUE)
plot_ethnic_roc(bellenguez_adsp_0318_pip01, col=col_roc_anno, title='max snp PIP > 0.1', plot=TRUE)
plot_ethnic_roc(bellenguez_adsp_0318_pip025, col=col_roc_anno, title='max snp PIP > 0.25', plot=TRUE)
