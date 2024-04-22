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
library(gridExtra)

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
  if("PRS_01" %in% colnames(df))
  {
    print("change PRS column name")
    #colnames(df) = gsub("_0\\.", "_", colnames(df))
    colnames(df) = gsub("PRS_", "P ",colnames(df))
    colnames(df) = gsub("P 0", "P 0.0",colnames(df))
    names(df)[names(df) == "P 1"] = 'P 0.1'
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


## new_anno_0318_24 with pip thres
bellenguez_adsp_0318_no_thres <- pre_process('/gpfs/commons/home/tlin/output/prs/new_anno_0318_24/bellenguez_adsp_reference/thres/prs_no_pip_thres.tsv')
bellenguez_adsp_0318_pip01 <- pre_process('/gpfs/commons/home/tlin/output/prs/new_anno_0318_24/bellenguez_adsp_reference/thres/prs_pip0.1.tsv')
bellenguez_adsp_0318_pip025 <- pre_process('/gpfs/commons/home/tlin/output/prs/new_anno_0318_24/bellenguez_adsp_reference/thres/prs_pip0.25.tsv')


pip_thres_0.1 <- pre_process('/gpfs/commons/home/tlin/output/prs/new_anno_0318_24/bellenguez_adsp_reference/thres/prs_pip_0.1.tsv')
pip_thres_0.2 <- pre_process('/gpfs/commons/home/tlin/output/prs/new_anno_0318_24/bellenguez_adsp_reference/thres/prs_pip_0.2.tsv')
pip_thres_0.25 <- pre_process('/gpfs/commons/home/tlin/output/prs/new_anno_0318_24/bellenguez_adsp_reference/thres/prs_pip_0.25.tsv')
pip_thres_0.3 <- pre_process('/gpfs/commons/home/tlin/output/prs/new_anno_0318_24/bellenguez_adsp_reference/thres/prs_pip_0.3.tsv')
pip_thres_0.4 <- pre_process('/gpfs/commons/home/tlin/output/prs/new_anno_0318_24/bellenguez_adsp_reference/thres/prs_pip_0.4.tsv')
pip_thres_0.5 <- pre_process('/gpfs/commons/home/tlin/output/prs/new_anno_0318_24/bellenguez_adsp_reference/thres/prs_pip_0.5.tsv')
pip_thres_0.6 <- pre_process('/gpfs/commons/home/tlin/output/prs/new_anno_0318_24/bellenguez_adsp_reference/thres/prs_pip_0.6.tsv')
pip_thres_0.7 <- pre_process('/gpfs/commons/home/tlin/output/prs/new_anno_0318_24/bellenguez_adsp_reference/thres/prs_pip_0.7.tsv')
pip_thres_0.8 <- pre_process('/gpfs/commons/home/tlin/output/prs/new_anno_0318_24/bellenguez_adsp_reference/thres/prs_pip_0.8.tsv')



pip_thres_susie <- pre_process('/gpfs/commons/home/tlin/output/prs/new_anno_0318_24/bellenguez_adsp_reference/thres/susie.tsv')
pip_thres_bl <- pre_process('/gpfs/commons/home/tlin/output/prs/new_anno_0318_24/bellenguez_adsp_reference/thres/baseline.tsv')
pip_thres_omics <- pre_process('/gpfs/commons/home/tlin/output/prs/new_anno_0318_24/bellenguez_adsp_reference/thres/omics.tsv')
pip_thres_omics_dl <- pre_process('/gpfs/commons/home/tlin/output/prs/new_anno_0318_24/bellenguez_adsp_reference/thres/omics_dl.tsv')

## new_anno_0318_14 with p thres
bellenguez_adsp_0318_susie <- pre_process('/gpfs/commons/home/tlin/output/prs/new_anno_0318_24/bellenguez_adsp_reference/p_thres/susie.tsv')
bellenguez_adsp_0318_bl <- pre_process('/gpfs/commons/home/tlin/output/prs/new_anno_0318_24/bellenguez_adsp_reference/p_thres/baseline.tsv')
bellenguez_adsp_0318_omics <- pre_process('/gpfs/commons/home/tlin/output/prs/new_anno_0318_24/bellenguez_adsp_reference/p_thres/omics.tsv')
bellenguez_adsp_0318_omics_dl <- pre_process('/gpfs/commons/home/tlin/output/prs/new_anno_0318_24/bellenguez_adsp_reference/p_thres/omics_dl.tsv')


## pT_CLUMPING
bellenguez_pT <- pre_process('/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/PRS/36k_ibd/bellenguez_pT_36k_ibd.tsv')
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
col_roc_polyfun_pthres <- list("P e5","P e4","P 0.001","P 0.01","P 0.1")
col_roc_pthres <- list("PRS_1","PRS_2","PRS_3","PRS_4","PRS_5","PRS_6","PRS_7","PRS_8")
#col_roc_polyfun_pthres <- list("PRS_e5","PRS_e4","PRS_001","PRS_01","PRS_1")

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




plot_ethnic_roc(pip_thres_susie, col=col_roc_pthres, title='SUSIE', plot=TRUE)
plot_ethnic_roc(pip_thres_bl, col=col_roc_pthres, title='baseline', plot=TRUE)
plot_ethnic_roc(pip_thres_omics, col=col_roc_pthres, title='omics', plot=TRUE)
plot_ethnic_roc(pip_thres_omics_dl, col=col_roc_pthres, title='omics_dl', plot=TRUE)

## test adsp reference panel without PIP thres, but with p thres

plot_ethnic_roc(bellenguez_adsp_0318_susie, col=col_roc_polyfun_pthres, title='susie', plot=TRUE)
plot_ethnic_roc(bellenguez_adsp_0318_bl, col=col_roc_polyfun_pthres, title='bl', plot=TRUE)
plot_ethnic_roc(bellenguez_adsp_0318_omics, col=col_roc_polyfun_pthres, title='omics', plot=TRUE)
plot_ethnic_roc(bellenguez_adsp_0318_omics_dl, col=col_roc_polyfun_pthres, title='omics_dl', plot=TRUE)


plot_ethnic_roc(bellenguez_pT, col=col_roc_polyfun_pthres, title='pT', plot=TRUE)



plot_ethnic_roc_facet(bellenguez_adsp_0318_susie ,bellenguez_adsp_0318_omics_dl, QC3 =bellenguez_pT,QC1name = "SuSiE", 
                      QC2name = "omics_dl",QC3name='pT',col = list("P e5","P 0.001","P 0.01","P 0.1"),boot_num = 50, title='bellenguez ',legendname = '')




### test corrlation


cor(bellenguez_adsp_0318_omics_dl$`P 0.1`,bellenguez_adsp_0318_omics_dl$Diagnosis )



print(cor(bellenguez_adsp_0318_omics_dl[, c('Diagnosis','P 0.1','P 0.01','P e5')]))


print(cor(bellenguez_adsp_0318_susie[, c('Diagnosis','P 0.1','P 0.01','P e5')]))


print(cor(bellenguez_adsp_0318_susie[bellenguez_adsp_0318_susie$predicted_ancestry =='EUR', c('Diagnosis','P 0.1','P 0.01','P e5')]))
print(cor(bellenguez_adsp_0318_susie[bellenguez_adsp_0318_susie$predicted_ancestry =='AMR', c('Diagnosis','P 0.1','P 0.01','P e5')]))
print(cor(bellenguez_adsp_0318_susie[bellenguez_adsp_0318_susie$predicted_ancestry =='AFR', c('Diagnosis','P 0.1','P 0.01','P e5')]))

roc(bellenguez_adsp_0318_susie[["Diagnosis"]],bellenguez_adsp_0318_susie$`P 0.1` , quiet=T)$auc

print(cor(bellenguez_pT[, c('Diagnosis','P 0.1','P 0.01','P e5')]))


bellenguez_adsp_0318_susie$Diagnosis <- factor(bellenguez_adsp_0318_susie$Diagnosis)

density_plot <- function(df,title,column)
{
  
  #df$Diagnosis <- factor(df$Diagnosis)
  df$Diagnosis <- factor(ifelse(df$Diagnosis == 0, "Control", "Case"))
  print(column)
  plot <- ggplot(df, aes_string(x=column, colour = 'Diagnosis', label = 'Diagnosis')) +
    geom_density() + theme_bw()  +labs(title = title, x = 'PRS',y= 'count') 
  return (plot)
}

pip_thres_densityplot <- function(df, title){
  plot1 <- density_plot(df, 'susie','PRS_susie')
  plot2 <- density_plot(df, 'bl','PRS_bl')
  plot3 <- density_plot(df, 'omics','PRS_bl_omics')
  plot4 <- density_plot(df, 'omics + dl','PRS_bl_omics_dl')
  grid.arrange(plot1, plot2, plot3,plot4, nrow = 2,top = textGrob(title,gp=gpar(fontsize=15,font=1)))
  
}

pip_thres_densityplot(pip_thres_0.1, 'pip_credibleset thres 0.1')
pip_thres_densityplot(extract_race(pip_thres_0.1,'EUR'), 'pip_credibleset thres 0.1 EUR')
pip_thres_densityplot(extract_race(pip_thres_0.1,'AMR'), 'pip_credibleset thres 0.1 AMR')
pip_thres_densityplot(extract_race(pip_thres_0.1,'AFR'), 'pip_credibleset thres 0.1 AFR')


pip_thres_densityplot(pip_thres_0.25, 'pip_credibleset thres 0.25')
pip_thres_densityplot(extract_race(pip_thres_0.25,'EUR'), 'pip_credibleset thres 0.25 EUR')
pip_thres_densityplot(extract_race(pip_thres_0.25,'AMR'), 'pip_credibleset thres 0.25 AMR')
pip_thres_densityplot(extract_race(pip_thres_0.25,'AFR'), 'pip_credibleset thres 0.25 AFR')



pip_thres_densityplot(pip_thres_0.8, 'pip_credibleset thres 0.8')
pip_thres_densityplot(extract_race(pip_thres_0.8,'EUR'), 'pip_credibleset thres 0.8 EUR')
pip_thres_densityplot(extract_race(pip_thres_0.8,'AMR'), 'pip_credibleset thres 0.8 AMR')
pip_thres_densityplot(extract_race(pip_thres_0.8,'AFR'), 'pip_credibleset thres 0.8 AFR')

pip_thres_densityplot(bellenguez_adsp_0318_pip01,"PIP thres 0.1 (all population)")
pip_thres_densityplot(extract_race(bellenguez_adsp_0318_pip01,'EUR'),"PIP thres 0.1 (EUR)")
pip_thres_densityplot(extract_race(bellenguez_adsp_0318_pip01,'AMR'),"PIP thres 0.1 (AMR)")
pip_thres_densityplot(extract_race(bellenguez_adsp_0318_pip01,'AFR'),"PIP thres 0.1 (AFR)")


pip_thres_densityplot(bellenguez_adsp_0318_no_thres,"no PIP thres (all population)")
pip_thres_densityplot(extract_race(bellenguez_adsp_0318_no_thres,'EUR'),"no PIP thres (EUR)")
pip_thres_densityplot(extract_race(bellenguez_adsp_0318_no_thres,'AMR'),"no PIP thres (AMR)")
pip_thres_densityplot(extract_race(bellenguez_adsp_0318_no_thres,'AFR'),"no PIP thres (AFR)")



pip_thres_densityplot(allele_flip, "allele_flip")

density_plot(plink_chr19, 'plink bl chr19','PRS')


density_plot(bellenguez_new_anno_0824_susie ,'','PRS')


density_plot(kunkle_adsp_no_apoe_qc ,'','PRS_5')
density_plot(kunkle_polyfun_pT ,'','PRS_5')



test_chr19 <- read.csv('/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_QC/annotated_hg37_plink_ibd/plink_prs/clump_bellenguez/test_chr19.prs.tsv', sep = '\t')
density_plot(test_chr19,'pt chr19', 'PRS')

hist(height$Height, height$PRS_05)


summary(height$Height)


ggplot(height[height$Height> 170.1,]$PRS_5) +boxplot()

subset_data1 <- subset(height, Height > 171)$PRS_5
subset_data2 <- subset(height, Height < 169)$PRS_5



# Create data frames for plotting
data1 <- data.frame(Height = rep("Height > 170.1", length(subset_data1)), PRS_5 = subset_data1)
data2 <- data.frame(Height = rep("Height < 169", length(subset_data2)), PRS_5 = subset_data2)
combined_data <- rbind(data1, data2)

# Create boxplot
ggplot(combined_data, aes(x = Height, y = PRS_5)) +
  geom_boxplot() +
  labs(x = "Height", y = "PRS_5") +
  theme_minimal()




## test wightman

wightman_pT_check<- read.csv('/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_QC/annotated_hg37_plink_ibd/plink_prs/clump_wightman/pT0.5_all_chrom.tsv', sep = '\t')

density_plot(wightman_pT_check,'whole_genome', 'SCORE_wholegenome')
density_plot(wightman_pT_check,'chr19', 'SCORE_19')

plots <- list()
for (chr in seq(1,22)) {
  plots[[chr]] <-density_plot(wightman_pT_check,sprintf("chr%d", chr), sprintf("SCORE_%d", chr))

}

grid.arrange(grobs = plots, nrow = 5,top = textGrob("in each chr", gp = gpar(fontsize = 15, font = 1)))



allele_flip <- pre_process('/gpfs/commons/home/tlin/output/prs/new_anno_0318_24/bellenguez_adsp_reference/allele_flip/allele_flip.prs.tsv')

cor(allele_flip$Diagnosis, allele_flip$PRS_bl_omics_dl)
plink_chr19<- pre_process('/gpfs/commons/home/tlin/output/prs/new_anno_0318_24/bellenguez_adsp_reference/thres/plink_output/sum_baseline_chr19.profilepheno.tsv')
plink_chr19= plink_chr19[plink_chr19$Age >= 65,]
density_plot(plink_chr19, 'plink bl chr19','PRS')
cor(plink_chr19$PRS, plink_chr19$Diagnosis)

cor(polypred_check$PRS, polypred_check$Diagnosis)



polypred_check <- read.csv('/gpfs/commons/home/tlin/output/prs/new_anno_0318_24/bellenguez_adsp_reference/check_result_bl_19.tsv.prspheno.tsv', sep ='\t')
polypred_check = polypred_check[polypred_check$Age >= 65,]
density_plot(polypred_check, 'polypred bl chr19','PRS')


cor(polypred_check$PRS, plink_chr19$PRS)
