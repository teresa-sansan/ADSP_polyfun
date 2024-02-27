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
## new plink 
kunkle_adsp_qc <- pre_process('/gpfs/commons/home/tlin/output/prs/new_plink/kunkle/kunkle_ADSP_qc_check.tsv')
kunkle_UKBB_qc <- pre_process('/gpfs/commons/home/tlin/output/prs/new_plink/kunkle/kunkle_ADSP_UKBB_qc.tsv')
kunkle_adsp_no_apoe_qc <- pre_process('/gpfs/commons/home/tlin/output/prs/new_plink/kunkle/kunkle_no_apoe_qc_check.tsv')

bellenguez_adsp_qc <- pre_process('/gpfs/commons/home/tlin/output/prs/new_plink/bellenguez/bellenguez_ADSP_qc_all.tsv')
bellenguez_UKBB_qc <- pre_process('/gpfs/commons/home/tlin/output/prs/new_plink/bellenguez/bellenguez_ADSP_UKBB_qc.tsv')
new_bellenguez_adsp <- pre_process('/gpfs/commons/home/tlin/output/prs/new_plink/bellenguez/new_sep22/ADSP.tsv')

wightman_UKBB_qc <- pre_process('/gpfs/commons/home/tlin/output/prs/new_plink/wightman/fixed_rsid1002/ADSP_qc_all.tsv')
jansen_adsp_qc <- pre_process('/gpfs/commons/home/tlin/output/prs/new_plink/jansen/ADSP_qc_all.tsv')

###  polyfun-Pred ----
kunkle_bl <- pre_process('/gpfs/commons/home/tlin/output/prs/polypred/kunkle/new_anno/bl.tsv')
kunkle_no_ml <- pre_process('/gpfs/commons/home/tlin/output/prs/polypred/kunkle/new_anno/no_ml.tsv')
kunkle_enformer <- pre_process('/gpfs/commons/home/tlin/output/prs/polypred/kunkle/new_anno/all_enformer.tsv')
kunkle_all_anno <- pre_process('/gpfs/commons/home/tlin/output/prs/polypred/kunkle/new_anno/all_anno.tsv')

bellenguez_bl <- pre_process('/gpfs/commons/home/tlin/output/prs/polypred/bellenguez/new_anno/bl.tsv')
bellenguez_no_ml <- pre_process('/gpfs/commons/home/tlin/output/prs/polypred/bellenguez/new_anno/no_ml.tsv')
bellenguez_enformer <- pre_process('/gpfs/commons/home/tlin/output/prs/polypred/bellenguez/new_anno/all_enformer.tsv')
bellenguez_all_anno <- pre_process('/gpfs/commons/home/tlin/output/prs/polypred/bellenguez/new_anno/all_anno.tsv')

wightman_susie <- pre_process('/gpfs/commons/home/tlin/output/prs/polypred/wightman/check_1003_susie.prs.tsv')
wightman_polypred <- pre_process('/gpfs/commons/home/tlin/output/prs/polypred/wightman/check_1003.prs.tsv')
wightman_bl <- pre_process('/gpfs/commons/home/tlin/output/prs/polypred/wightman/check_1003_bl.prs.tsv')
wightman_fix_convergence <- pre_process('/gpfs/commons/home/tlin/output/prs/polypred/wightman/check_1003_fixed_convergence.tsv')

wightman_no_ml <- pre_process('/gpfs/commons/home/tlin/output/prs/polypred/wightman/no_ml.tsv')
wightman_all_anno <- pre_process('/gpfs/commons/home/tlin/output/prs/polypred/wightman/update_all+enformer.tsv')
wightman_enformer <- pre_process('/gpfs/commons/home/tlin/output/prs/polypred/wightman/enformer.tsv')
wightman_all_except_enformer <- pre_process('/gpfs/commons/home/tlin/output/prs/polypred/wightman/all_except_enformer.tsv')
wightman_glasslab <- pre_process('/gpfs/commons/home/tlin/output/prs/polypred/wightman/glasslab.tsv')


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


## Extract specific race ----- 
extract_race <-function (df,race){
  one_race = df[df$final_population == race,]
  return(one_race)
}

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

## process PRS name ------
process_prs_col_name <- function(df) {
  if ('PRS10' %in% df$PRS) {
    df$PRS <- sub('^PRS', 'max snp', df$PRS)
    if ('PRS3' %in% df$PRS) {
      df$PRS <- factor(df$PRS, level = c('max snp 10', 'max snp 7', 'max snp 5', 'max snp 3', 'max snp 1'))
    } else if ('PRS2' %in% df$PRS) {
      df$PRS <- factor(df$PRS, level = c('max snp 5', 'max snp 2', 'max snp 1'))
    } else {
      df$PRS <- factor(df$PRS, level = c('max snp 10', 'max snp 5', 'max snp 1'))
    }
  } else if ('PRS_5' %in% df$PRS) {
    df$PRS <- sub('^PRS_', 'p = 0.', df$PRS)
    df[df$PRS == 'p = 0.e5', 'PRS'] <- 'p = 1e-5'
    df$PRS <- factor(df$PRS, level = c('p = 1e-5', 'p = 0.001', 'p = 0.005', 'p = 0.01', 'p = 0.05', 'p = 0.1', 'p = 0.5'))
  }
  return(df)
}

## AUC ------
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
## If plot = T(default), will return a plot with AUC using different method/threshold, if false, will return a list for other purpose. 

roc_result <-function(df, column_for_roc = col_roc_E5){
  auc_list=list()
  auc_df = data.frame("PRS" = unlist(column_for_roc), "auc" = 0)
  for (i in 1:length(column_for_roc)){
    col = column_for_roc[[i]]
    roc_formula = roc(df[["Diagnosis"]],df[[col]], quiet=T)
    auc_value <- roc(df[["Diagnosis"]],df[[col]], quiet=T)$auc
    auc_list[i] = as.numeric(auc_value)
  }
  auc_df$auc = unlist(auc_list)
  return(auc_df)
}
roc_result_boot <-function(df,column_for_roc=col_roc_E5, boot_num=50, mean=FALSE){
  CI_roc = data.frame(matrix(ncol = 4, nrow = 0))
  colnames(CI_roc) = c("PRS","boot_CI_lower","boot_mean","boot_CI_upper")
  for (i in 1:length(column_for_roc)){
    col = column_for_roc[[i]]
    roc_formula <- roc(df[["Diagnosis"]],df[[col]], quiet=T)
    set.seed(10)
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
## remember to change x lim
plot_ethnic_roc <- function(df, col=col_roc_E5,boot_num=50, boot=TRUE, plot=FALSE, title=''){
  output_df <- data.frame(matrix(ncol = 0, nrow = length(col)*3))
  output_df$PRS =  rep(unlist(col),1)
  output_df$ethnicity = c(rep("EUR",length(col)),rep("AFR",length(col)),rep("AMR",length(col)))
  output_df$ethnicity <- factor(output_df$ethnicity,levels = c("EUR","AFR","AMR"))
  
  if (boot != TRUE) {
    races <- c('EUR', 'AFR', 'AMR')
    output_df$auc <- unlist(lapply(races, function(race) roc_result(extract_race(df, race), column_for_roc = col)$auc))
  } else {
    print(paste("boot time =", boot_num))
    races <- c('EUR', 'AFR', 'AMR')
    boot_results <- lapply(races, function(race) roc_result_boot(extract_race(df, race), column_for_roc = col, boot_num = boot_num))
    
    output_df <- do.call(rbind, lapply(boot_results, function(result) c(result$boot_mean, result$boot_CI_lower, result$boot_CI_upper)))
    output_df[, c("auc", "boot_CI_lower", "boot_CI_upper")] <- as.numeric(output_df[, c("auc", "boot_CI_lower", "boot_CI_upper")])
  }
  
  if(plot != FALSE){
    output_df = process_prs_col_name(output_df)
    plot <- ggplot(data = output_df, aes(x=auc, y = PRS))+ geom_point(size=3,alpha=0.9,position = position_dodge(width = 0.7), color='darkblue')+
      facet_wrap(~ethnicity, ncol=1)+ xlab('AUC')+ ggtitle(title)+xlim(0.3, 0.72)+theme_bw()
    if (boot){  
      plot <- plot + geom_errorbar(aes(xmin=boot_CI_lower, xmax=boot_CI_upper),position=position_dodge(width=0.7), width=.1,alpha=0.5,
                                   color='darkblue',show.legend = FALSE) 
    }## plot error bar or not
    print(plot)
  }
  return(output_df)
}

## fixed legend
## The boolean APOE operator is to see whether you want to plot the effect of only using 5(6, depending on QC or not) APOE allele to make prediction. 
plot_ethnic_roc_facet <- function(QC1, QC2, QC3,QC4=data.frame(),QC5=data.frame(),QC6=data.frame(), col=col_roc_E5, title=' ',
                                  QC1name="", QC2name="",QC3name="",QC4name=' ',QC5name=' ',QC6name=' ',boot=TRUE, 
                                  boot_num=50, legendname=FALSE, APOE=FALSE){
  dfs <- list(QC1, QC2, QC3, QC4, QC5, QC6)
  dfsname <- list(QC1name, QC2name, QC3name, QC4name, QC5name, QC6name)
  
  roc_add_anno <- function(df, col, boot, bootnum, anno_name){
    df = plot_ethnic_roc(df, col=col,boot=boot,boot_num=boot_num)
    df$qc_status = anno_name
    return(df)
  } 
  result <- do.call(rbind, mapply(function(df, name) {
    if (nrow(df) != 0) {
      name <- factor(name, levels = unique(dfsname))
      df_new <- roc_add_anno(df, col, boot, boot_num, name)
      return(df_new)
    }}, dfs, dfsname, SIMPLIFY = FALSE))
  
  result = process_prs_col_name(result)
  print(result)
  
  result$PRS <- str_replace(result$PRS,'PRS','max snp ')
  result$PRS <- factor(result$PRS, levels=unique(result$PRS))

  plot <- ggplot(data = result, aes(x=auc, y = PRS, color = qc_status,))+
    geom_point(size=3,alpha=0.7,position = position_dodge(width = 0.7))+
    facet_wrap(~ethnicity, ncol=1)+
    xlab('AUC')+ ggtitle(title)+xlim(0.45, 0.7)+
    theme_bw()
  
  if (boot){  
    plot <- plot + geom_errorbar(aes(xmin=boot_CI_lower, xmax=boot_CI_upper),position=position_dodge(width=0.7), width=.1,alpha=0.5,show.legend = FALSE) 
  }## plot error bar or not
  if(legendname != FALSE){ 
    plot = plot +guides(col=guide_legend(legendname))
  } ## change legend name or not (depends on usage)
  if(class(APOE) != "logical"){
    APOE_ci = plot_ethnic_roc(APOE, col="PRS", plot='F',boot=boot, boot_num=boot_num)
    plot = plot +  geom_vline(data = APOE_ci, aes(xintercept = auc,color='only APOE'),color="black",lwd=0.5, alpha=0.7,linetype="dashed")  ## move color into aes will generate legend automatically.
  } ## if we want to plot APOE (line) ## should do bootstrap automatically if it is specified 
  return(plot)
}

plot_auc_facet_all_sumstat <- function(a,b,c,d,legendname = FALSE){
  prow <- plot_grid(a+ theme(legend.position="none"), 
                    b+ theme(legend.position="none",axis.text.y = element_blank())+ ylab(NULL), 
                    c+ theme(legend.position="none",axis.text.y = element_blank())+ylab(NULL),
                    d+ theme(legend.position="none",axis.text.y = element_blank())+ylab(NULL),
                    ncol = 3, nrow = 1)
  legend <- get_legend(
    d +guides(col=guide_legend(legendname, nrow=1))+
      theme(legend.position = "bottom")
  )
  plot_grid(prow, legend, ncol=1,rel_heights=c(3,.4))
  
}

## R2 functions ----
## calculate pseudo rsquare 
rsq_formula <- function(formula, data, indices=FALSE) {    
  d <- data[indices,] # allows boot to select sample
  fit <- glm(formula = formula,family = 'binomial', data = d)
  return(RsqGLM(fit, plot = FALSE)$Nagelkerke)
}## this is for bootstrapping

## relative R2 w. bootstrap -----
##resample
log_reg <- function(df,prs,main_title, plot = TRUE, legend="PRS_", boot_num = FALSE, replace="pT = 0."){
  set.seed(50)
  mod1 <- glm(Diagnosis ~ Sex + Age, data= df, family=binomial) ## first create a formula only using covariates
  mod1_R2 <- RsqGLM(mod1, plot=FALSE)$Nagelkerke
  for (i in 1:length(prs)){
    frm <- as.formula(paste("Diagnosis ~ Sex + Age+", prs[[i]])) ## add the PRS you wanted
    if(boot_num != FALSE){              
      mod <- boot(data=df, statistic=rsq_formula,
                  R=boot_num, formula=frm)
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
    plot(log_output$boot_mean, ylim  = c(min(log_output$null_R2[1], min(log_output$boot_CI_lower))*0.1, max(log_output$boot_mean[1], max(log_output$boot_CT_upper)*1.5)),
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

## if not doing bootstraping here, we are calculating partial R2
plot_ethnic_R2 <- function(df, col, boot_num, title=' ',replace='', plot= FALSE){
  if(class(boot_num) == 'numeric'){
    print(paste("running EUR with bootstrapping ", boot_num, ' times'))
    EUR = log_reg(extract_race(df,'EUR'), col,  paste(title, ", EUR"), plot = plot, legend="PRS",boot_num = boot_num, replace=replace)
    print(paste("running AFR with bootstrapping ", boot_num, ' times'))
    AFR = log_reg(extract_race(df,'AFR'), col, paste(title, ", AFR"), plot = plot, legend="PRS",boot_num = boot_num, replace=replace)
    print(paste("running AMR with bootstrapping ", boot_num, ' times'))
    AMR = log_reg(extract_race(df,'AMR'), col,  paste(title, ", AMR"), plot = plot, legend="PRS",boot_num = boot_num, replace=replace)
  } ## do bootstrapping
  else{
    EUR = log_reg_partial_cov(extract_race(df,'EUR'), col)
    AFR = log_reg_partial_cov(extract_race(df,'AFR'), col)
    AMR = log_reg_partial_cov(extract_race(df,'AMR'), col)
     ## calculate partial R2
  }
  EUR$ethnicity="EUR"
  AFR$ethnicity="AFR"
  AMR$ethnicity="AMR"
  df = rbind(EUR, AFR, AMR)
  df$ethnicity <- factor(df$ethnicity, levels = c("EUR","AFR","AMR"))
  if(plot == FALSE){
    return(df)
  }
}

## plot facet 
## plot ethnic facut plot for one sumstat

plot_ethnic_R2_facet <- function(QC1, QC2, QC3, col=col_roc_E5, boot_num=50, title=' ',QC1name="QC_on_base", QC2name="QC_on_base+variants",QC3name="no_qc", legendname=FALSE){
  QC1 = plot_ethnic_R2(QC1, col, boot_num)
  QC2 = plot_ethnic_R2(QC2, col, boot_num)
  QC1$qc_status = QC1name
  QC2$qc_status = QC2name
  if(class(QC3) != 'logical'){
    QC3 = plot_ethnic_R2(QC3, col, boot_num)
    QC3$qc_status = QC3name
    all = rbind(QC1,QC2,QC3)
  }else
    all = rbind(QC1, QC2)
  all = process_prs_col_name(all)
  add_xlegend = '' 
  if("partial_R2" %in% colnames(all)){
    print('R2 here is partial R2')
    all$boot_mean = all$partial_R2
    all$boot_CI_lower = 0
    all$boot_CI_upper = 0
    add_xlegend = 'Partial ' 
  }
  all$boot_mean = all$boot_mean * 100
  all$boot_CI_upper = all$boot_CI_upper * 100
  all$boot_CI_lower = all$boot_CI_lower * 100
  plot <- ggplot(data = all, aes(x= boot_mean, y = PRS, color = qc_status))+
    geom_point(size=2, alpha=0.7,position = position_dodge(width = 0.7))+
    geom_pointrange(aes(xmin=boot_CI_lower, xmax=boot_CI_upper), linetype="dotted",position=position_dodge(width=0.7),show.legend = FALSE) +
    facet_wrap(~ethnicity, ncol=1)+
    xlab(paste(add_xlegend,'R squared (%)'))+ ggtitle(title)+xlim(min(all$boot_CI_lower), max(all$boot_CI_upper)) +
    theme_bw()
  if(legendname != FALSE){  ## testout new legend title
    plot = plot +guides(col=guide_legend(legendname))
  }
  print(all)
  return(plot)
}

plot_ethnic_R2_facet_race <- function(QC1, QC2, QC3, col=col_roc_polypred3, boot_num=FALSE, title=' ',QC1name="bl", QC2name="bl+enformer",QC3name="all", legendname=FALSE){
  QC1 = plot_ethnic_R2(QC1, col, boot_num)
  QC2 = plot_ethnic_R2(QC2, col, boot_num)
  QC1$qc_status = QC1name
  QC2$qc_status = QC2name
  if(class(QC3) != 'logical'){
    QC3 = plot_ethnic_R2(QC3, col, boot_num)
    QC3$qc_status = QC3name
    all = rbind(QC1,QC2,QC3)
  }else
    all = rbind(QC1, QC2)
  all=all[all$PRS!="PRS5",]
  all = process_prs_col_name(all)
  
  add_xlegend = '' 
  if("partial_R2" %in% colnames(all)){
    print('R2 here is partial R2')
    all$sex_age = all$partial_R2
    all$sex_age_pc = all$partial_R2_pc
    all$sex_pc = all$partial_R2_noage
    add_xlegend = 'Partial ' 
  }
  plot <- ggplot(data = all, aes(x= sex_age_pc, y = PRS, color = qc_status))+
    geom_point(size=2, alpha=0.7,position = position_dodge(width = 0.7))+
    facet_wrap(~ethnicity, ncol=1)+guides(color = guide_legend(title = "Annotations"))+
    xlab(paste(add_xlegend,'R squared (%), + 5PCs'))+ ggtitle(title)+
    theme_bw()
  if(legendname != FALSE){  ## testout new legend title
    plot = plot +guides(col=guide_legend(legendname))
  }
  print(all)
  return(plot)
}

## plot ethnic facut plot for multiple sumstat
plot_R2_facet_allsumstat<- function(s1_qc1, s1_qc2, s1_qc3,s2_qc1, s2_qc2, s2_qc3,s3_qc1, s3_qc2, s3_qc3,col = col_roc_E5, QC1name="QC_on_base", QC2name="QC_on_base+variants",QC3name="no_qc",boot_num=50, legendname=legendname){
  a = plot_ethnic_R2_facet(s1_qc1, s1_qc2, QC3=s1_qc3, col=col, boot_num=boot_num, 'Kunkle', QC1name, QC2name,QC3name,legendname=legendname)
  b = plot_ethnic_R2_facet(s2_qc1, s2_qc2, QC3=s2_qc3, col=col, boot_num=boot_num, 'Bellenguez',QC1name, QC2name,QC3name,legendname=legendname)
  c = plot_ethnic_R2_facet(s3_qc1, s3_qc2, QC3=s3_qc3, col=col, boot_num=boot_num,'Wightman', QC1name, QC2name,QC3name,legendname=legendname)

  prow <- plot_grid(a+ theme(legend.position="none"), 
                    b+ theme(legend.position="none",axis.text.y = element_blank())+ ylab(NULL), 
                    c+ theme(legend.position="none",axis.text.y = element_blank())+ylab(NULL),
                    ncol = 3, nrow = 1)
  legend_b <- get_legend(
    c+ guides(color = guide_legend(legendname, nrow = 1)) +
      theme(legend.position = "bottom")
  )
  print(prow)
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


## partial R2 -----
## check if I can fixed the error in 'numeric 'envir' arg not of length one ' when using not specified df in rsq function.
## code is from https://github.com/cran/rsq/blob/master/R/rsq.R
## directly apply function will be faster
## Gave up on bootstrap for partial R2

log_reg_partial <- function(df,col){
  modR <- glm(Diagnosis ~ Sex + Age, data= df, family=binomial) ## first create a reduced model only using covariates
  modR_R2 <- RsqGLM(modR, plot=FALSE)$Nagelkerke
  partial_R2 = rep(NA, length(col))
  for (i in 1:length(col)){
    frm <- as.formula(paste("Diagnosis ~ Sex + Age + ", col[[i]])) ## add the PRS you wanted
    modF <- glm(formula = frm,family = 'binomial', data = df) ## with PRS
    modF_R2 <-  RsqGLM(modF, plot=FALSE)$Nagelkerke 
    partialR2 <-  1- ((1-modF_R2 ) / (1-modR_R2))
    partial_R2[i] = partialR2*100
  }
  PRS =  array(unlist(col))
  return(data.frame(PRS,partial_R2))
}


## need to re-write to make it more efficient
log_reg_partial_cov <- function(df,col){
  modR <- glm(Diagnosis ~ Sex + Age, data= df, family=binomial) ## first create a reduced model only using covariates
  modR_R2 <- RsqGLM(modR, plot=FALSE)$Nagelkerke
  
  modR_pc <- glm(Diagnosis ~ Sex + Age +X1+X2+X3+X4+X5, data= df, family=binomial) ## first create a reduced model only using covariates
  modR_R2_pc <- RsqGLM(modR_pc, plot=FALSE)$Nagelkerke
  
  modR_noage <- glm(Diagnosis ~ Sex +X1+X2+X3+X4+X5, data= df, family=binomial) ## first create a reduced model only using covariates
  modR_R2_noage <- RsqGLM(modR_noage, plot=FALSE)$Nagelkerke
  
  partial_R2 = rep(NA, length(col))
  partial_R2_pc = rep(NA, length(col))
  partial_R2_noage = rep(NA, length(col))
  
  for (i in 1:length(col)){
    frm <- as.formula(paste("Diagnosis ~ Sex + Age + ", col[[i]])) ## add the PRS you wanted
    modF <- glm(formula = frm,family = 'binomial', data = df) ## with PRS
    modF_R2 <-  RsqGLM(modF, plot=FALSE)$Nagelkerke 
    partialR2 <-  1- ((1-modF_R2 ) / (1-modR_R2))
    partial_R2[i] = partialR2*100
  
    frm <- as.formula(paste("Diagnosis ~ Sex + Age + X1 + X2 + X3 + X4 + X5 + ", col[[i]])) ## add the PRS you wanted
    modF <- glm(formula = frm,family = 'binomial', data = df) ## with PRS
    modF_R2 <-  RsqGLM(modF, plot=FALSE)$Nagelkerke 
    partialR2_pc <-  1- ((1-modF_R2 ) / (1-modR_R2_pc))
    partial_R2_pc[i] = partialR2_pc*100
  
    frm <- as.formula(paste("Diagnosis ~ Sex  + X1 + X2 + X3 + X4 + X5 + ", col[[i]])) ## add the PRS you wanted
    modF <- glm(formula = frm,family = 'binomial', data = df) ## with PRS
    modF_R2 <-  RsqGLM(modF, plot=FALSE)$Nagelkerke 
    partialR2_noage <-  1- ((1-modF_R2 ) / (1-modR_R2_noage))
    partial_R2_noage[i] = partialR2_noage*100
  }
  
  PRS =  array(unlist(col))
  return(data.frame(PRS,partial_R2,partial_R2_pc, partial_R2_noage))
}
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
plot_ethnic_roc_facet(kunkle_susie, kunkle_polypred, kunkle_bl, QC1name = "none",
                      QC2name = "all",QC3name='baseline', col = col_roc_polypred3,boot_num = 50, title='kunkle',legendname = 'annotation')

plot_ethnic_roc_facet(bellenguez_susie, bellenguez_polypred, bellenguez_bl, QC1name = "none",
                      QC2name = "all",QC3name='baseline', col = col_roc_polypred3,boot_num = 50, title='bellenguez', legendname = 'annotation')

plot_ethnic_R2_facet(wightman_susie_converged, wightman_bl_converged, wightman_polypred_converged, col = col_roc_polypred3, QC1name='none', QC2name = 'baseline', QC3name = 'all', title = 'wightman', legendname = 'annotations' )

plot_ethnic_roc_facet(wightman_polypred_converged, wightman_bl_converged, wightman_susie_converged, col = col_roc_polypred3, QC1name='all', QC2name = 'baseline', QC3name = 'none', title = 'wightman', legendname = 'annotations' )

plot_ethnic_roc_facet(wightman_susie, wightman_bl, wightman_polypred, col = col_roc_polypred3, QC1name='none', QC2name = 'baseline', QC3name = 'all', title = 'wightman', legendname = 'annotations' )

plot_ethnic_R2(wightman_glasslab, boot_num=FALSE, col=col_roc_polypred3)


## PRSCS ------
# Check the levels of the 'Diagnosis' variable
# Convert 'Diagnosis' from factor to numeric
polypred_PRSCS_NOT0$Diagnosis <- as.numeric(polypred_PRSCS_NOT0$Diagnosis)
polypred_prscs_noEnformer$Diagnosis <- as.numeric(polypred_prscs_noEnformer$Diagnosis)

plot_ethnic_roc_facet(polypred_PRSCS_NOT0, polyfun_PRSCS_NOT0, wightman_UKBB_qc,title='wightman',
                      QC1name = "polyfun_PRSCS(PIP >0)", QC2name = "polyfun_sumstat_PRSCS(PIP > 0)", QC3name='genomewide_plink',QC4name='genomewide_plink',
                      legendname = "method", boot_num=50)

plot_ethnic_roc_facet(PRSCS, polypred_PRSCS, polypred_plink,wightman_UKBB_qc,title='wightman',
                      QC1name = "PRSCS", QC2name = "polyfun_PRSCS", QC3name='polyfun_plink',QC4name='genomewide_plink',
                      legendname = "method", boot_num=2)

plot_ethnic_roc_facet(wightman_bl,wightman_enformer, wightman_all_except_enformer,wightman_all_anno, col = list("PRS5","PRS10"),
                      QC1name = "baseline", QC2name = "bl+enformer", QC3name='all but enformer',QC4name = 'all anno',
                      title = 'wightman', legendname = 'annotations' )

plot_ethnic_roc_facet(wightman_bl, wightman_enformer, wightman_all_anno, wightman_no_ml, col = list("PRS5","PRS10"),  
                      QC1name='bl',QC4name = 'no ml',QC2name = 'bl+enformer', QC3name = 'all anno',
                      title = 'wightman', legendname = 'annotations', boot_num = 50)

plot_ethnic_roc_facet(kunkle_bl, kunkle_no_ml, kunkle_enformer, kunkle_all_anno, col = list("PRS5","PRS2"),  
                      QC1name='bl',QC2name = 'no ml',QC3name = 'bl+enformer', QC4name = 'all anno',
                      title = 'kunkle', legendname = 'annotations', boot_num = 50)

plot_ethnic_roc_facet(bellenguez_bl, bellenguez_no_ml, bellenguez_enformer, bellenguez_all_anno, col =list("PRS5","PRS2"),  
                      QC1name='bl',QC2name = 'no ml',QC3name = 'bl+enformer', QC4name = 'all anno',
                      title = 'bellenguez', legendname = 'annotations' )

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
  xlab('AUC')+ ggtitle('Wightman')+xlim(0.3, 0.72)+theme_bw() +facet_wrap(~ethnicity, ncol=1)+
  geom_errorbar(aes(xmin=boot_CI_lower, xmax=boot_CI_upper),position=position_dodge(width=0.7), width=.1,alpha=0.5,color='darkblue',show.legend = FALSE) 


## partial r2
data_list <- list(kunkle_polypred, bellenguez_polypred, wightman_polypred, jansen_polypred)
lapply(data_list, apply_log_reg_partial)

plot_R2_facet_allsumstat(kunkle_polypred,  kunkle_polyfun_plink_no_cpT, FALSE,
                         bellenguez_adsp_qc, bellenguez_polyfun_plink_no_cpT,FALSE,
                         wightman_adsp_qc, wightman_polyfun_plink_no_cpT,FALSE,
                         QC1name = 'Summary Stat', QC2name = 'PolyFun',legendname = 'Effect size', boot_num=FALSE)

## bellenguez
#df,prs,main_title, plot = TRUE, legend="PRS_", boot_num = FALSE, replace="pT = 0."
plot_R2_cT(bellenguez_cT, bellenguez_qc, bellenguez_qc_on_base, bellenguez_qc_on_target, 'bellenguez_updateRSID', col_roc, plot=TRUE, boot_num = 50)
plot_R2_cT(bellenguez_cT0224, bellenguez_qc0224, bellenguez_qc_on_base0224, bellenguez_qc_on_target0224, 'bellenguez_fixed0224', col_roc_E5, boot_num=50 )

## interested SNP
mod1 <- glm(Diagnosis ~ Sex + Age + PRS, data= bellenguez_interested, family=binomial)
RsqGLM(mod1, plot=FALSE)$Nagelkerke

plot_R2_facet_allsumstat(kunkle_adsp_qc,  kunkle_susie, FALSE,
                         bellenguez_polypred, bellenguez_susie,FALSE,
                         wightman_polypred, wightman_polyfun_plink_no_cpT,FALSE,
                         QC1name = 'Summary Stat', QC2name = 'PolyFun',legendname = 'Effect size', boot_num=FALSE)
## new bellenguez 0922
plot_ethnic_roc_facet(bellenguez_adsp,new_bellenguez_adsp, FALSE, QC1name = 'old bellenugez',QC2name = 'bellenguez 0922', boot_num = 50, legendname = 'bellenguez version')

##bellenguez_interested_SNP 
df_list <- list(bellenguez_interested, bellenguez_interest_max,bellenguez_interest_min,
                extract_eur(bellenguez_interested),extract_eur(bellenguez_interest_max),extract_eur(bellenguez_interest_min),
                extract_afr(bellenguez_interested),extract_afr(bellenguez_interest_max),extract_afr(bellenguez_interest_min),
                extract_amr(bellenguez_interested),extract_amr(bellenguez_interest_max),extract_amr(bellenguez_interest_min))

for (i in df_list){
  print(auc(roc(Diagnosis~PRS, data = i, quiet=T))[1])
}

## APOE
plot_ethnic_roc(kunkle_withoutAPOE, "kunkle remove APOE region", col_roc_E5)
col_roc <- list("only2SNP","no2SNP","no2SNP_qc")

kunkle_APOE_roc <- roc(Diagnosis ~ only2SNP+no2SNP+no2SNP_qc, data = kunkle_APOE)
auc(roc(Diagnosis~only2SNP, data = kunkle_APOE)) ##0.621
auc(roc(Diagnosis~only2SNP, data = extract_eur(kunkle_APOE))) ##0.6187

for (i in list(extract_eur(kunkle_APOE),extract_afr(kunkle_APOE),extract_amr(kunkle_APOE))){
  print(auc(roc(Diagnosis~no2SNP, data = i, quiet=T))[1])
}


##kunkle

plot_ethnic_roc_facet(kunkle_susie, kunkle_bl, kunkle_polypred, col = col_roc_polypred3, title= 'kunkle',
                      QC1name = 'SuSiE', QC2name = 'PolyFun (BL)',QC3name='PolyFun (All anno)',
                      legendname = 'Annotation'
)

wightman_polypred$PRS = wightman_polypred$PRS10
plot_ethnic_R2_facet(wightman_susie_max10_polypred, wightman_polypred, FALSE, col = 'PRS', title= 'wightman',
                     QC1name = 'no annotations(SuSiE)', QC2name = 'Baseline annotations', QC3name = 'All annotations',
                     legendname ='Annotation',boot_num = FALSE)

bellenguez_polypred$PRS = bellenguez_polypred$PRS10
plot_ethnic_R2_facet(bellenguez_susie, bellenguez_polypred,FALSE,col ='PRS', title= 'bellenguez (old_plink)',
                      QC1name = 'SuSiE', QC2name = 'Polyfun',
                      legendname = 'Annotation',boot_num=FALSE
)
plot_ethnic_roc_facet(bellenguez_susie, bellenguez_polypred,FALSE,col =col_roc_polypred3, title= 'bellenguez (old_plink)',
                      QC1name = 'SuSiE', QC2name = 'Polyfun',
                      legendname = 'Annotation'
)

plot_ethnic_roc_facet(bellenguez_susie, bellenguez_polypred,FALSE,col ='PRS', title= 'bellenguez (old_plink)',
                      QC1name = 'SuSiE', QC2name = 'Polyfun',
                      legendname = 'Annotation',boot_num = FALSE
)


plot_ethnic_roc_facet(bellenguez_susie, bellenguez_polypred,FALSE,col =col_roc_polypred3, title= 'bellenguez (old_plink)',
                      QC1name = 'SuSiE', QC2name = 'Polyfun_pred',
                      legendname = 'Annotation'
)

plot_ethnic_roc_facet(bellenguez_polypred_new, bellenguez_susie, FALSE, boot_num = 50, col='PRS',
                      title='bellenguez', QC1name = 'PolyFun_Pred', QC2name = 'SuSiE', legendname = 'Tools')

plot_ethnic_roc_facet(wightman_susie_max10, wightman_polypred, FALSE, boot_num = 50, col='PRS',
                      title='wightman', QC1name = 'PolyFun_Pred', QC2name = 'SuSiE', legendname = 'Tools')



plot_ethnic_roc_facet(PRSCS_36k_ibd_wightman,inf_0214, inf_inter2_0219,  boot_num = 50, col='PRS',
                      title='wightman', QC1name = 'PRSCS', QC2name = 'inf', QC3name='inf_intersection2',legendname = 'Tools')


plot_ethnic_R2_facet(wightman_susie_max10, wightman_polypred,FALSE, col ='PRS', title= 'wightman',
                     QC1name = 'SuSiE', QC2name = 'Polyfun_pred',
                     legendname = 'Annotation',boot_num=FALSE
)


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

plot_R2_facet_allsumstat(kunkle_adsp_qc,  kunkle_polyfun_plink_no_cpT, FALSE,
                         bellenguez_adsp_qc, bellenguez_polyfun_plink_no_cpT,FALSE,
                         wightman_adsp_qc, wightman_polyfun_plink_no_cpT,FALSE,
                         QC1name = 'Summary Stat', QC2name = 'PolyFun',legendname = 'Effect size', boot_num=FALSE)

## other combination ------
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



for ( i in list(kunkle_jk, bellenguez_jk, wightman_jk, jansen_jk)){
  print(log_reg_partial(extract_eur(i), col_jk))
}

## kunkle
## remove all APOE
plot_ethnic_R2(kunkle_withoutAPOE, col_roc_E5,"Kunkle without APOE", 5)
plot_ethnic_R2(kunkle_withoutAPOE_qc, col_roc_E5,"Kunkle without APOE, QC on sumstat_variant", 50)
plot_ethnic_roc_facet(kunkle_cT,kunkle_withoutAPOE_qc,FALSE,title='kunkle',
                      QC1name= "with APOE", QC2name="without APOE, QC")

## remove two APOE allele
col_roc_APOE = list("only2SNP","no2SNP","no2SNP_qc")
kunkle_cT_log<-log_reg(kunkle_APOE, col_roc_APOE, "kunkle_cT", plot=FALSE)   
kunkle_cT_qc_log_eur<-log_reg(extract_eur(kunkle_APOE), col_roc_APOE, "kunkle QC", plot=FALSE)
kunkle_cT_qc_log_afr<-log_reg(extract_afr(kunkle_APOE), col_roc_APOE, "kunkle QC", plot=FALSE)
kunkle_cT_qc_log_amr<-log_reg(kunkle_APOE[kunkle_APOE$final_population == "AMR",], col_roc_APOE, "kunkle QC", plot=FALSE)

## wightman
plot_ethnic_R2(wightman_cT, col_roc_E5, 'Wightman', 50)

##susie 
plot_ethnic_R2(susie, col_roc_polypred, "bellenguez_susie",50)
plot_ethnic_R2(kunkle_susie, col_roc_polypred, "kunkle_susie",50)

RsqGLM(glm(Diagnosis~PRS10+Sex+Age,family = 'binomial', data = kunkle_susie), plot=FALSE)$Nagelkerke ## 0.14
RsqGLM(glm(Diagnosis~PRS1+Sex+Age,family = 'binomial', data = kunkle_susie), plot=FALSE)$Nagelkerke ##0.139



## test age threshold -----

control_age_roc <- function(df, age, col=col_roc_E5,boot_num=50, boot=TRUE, plot=FALSE, title=''){
  output_df <- data.frame(matrix(ncol = 0, nrow = length(col)*3))
  output_df$PRS =  rep(unlist(col),1)
  output_df$ethnicity = c(rep("EUR",length(col)),rep("AFR",length(col)),rep("AMR",length(col)))
  output_df$ethnicity <- factor(output_df$ethnicity,levels = c("EUR","AFR","AMR"))
  print('age segregation')
  EUR = roc_result_boot(segregate_control_by_age(extract_eur(df),age), column_for_roc = col, boot_num = boot_num)
  AFR = roc_result_boot(segregate_control_by_age(extract_afr(df),age), column_for_roc = col, boot_num = boot_num)
  AMR = roc_result_boot(segregate_control_by_age(extract_amr(df),age), column_for_roc = col, boot_num = boot_num)
  output_df$auc = append(append(unlist(EUR$boot_mean),unlist(AFR$boot_mean)),unlist(AMR$boot_mean))
  output_df$boot_CI_lower = append(append(unlist(EUR$boot_CI_lower),unlist(AFR$boot_CI_lower)),unlist(AMR$boot_CI_lower))
  output_df$boot_CI_upper = append(append(unlist(EUR$boot_CI_upper),unlist(AFR$boot_CI_upper)),unlist(AMR$boot_CI_upper))
  output_df[c("auc","boot_CI_lower","boot_CI_upper")] <- sapply(output_df[c("auc","boot_CI_lower","boot_CI_upper")],as.numeric)
  output_df$age = age
  
  if(plot != FALSE){
    output_df = process_prs_col_name(output_df)
    plot <- ggplot(data = output_df, aes(x=auc, y = PRS))+
      geom_point(size=3,alpha=0.9,position = position_dodge(width = 0.7), color='darkblue')+
      facet_wrap(~ethnicity, ncol=1)+
      xlab('AUC')+ ggtitle(paste(title, "age = ", as.character(age)))+xlim(0.45, max(output_df$boot_CI_upper)*1.1)+theme_bw()
    if (boot == TRUE){  
      plot <- plot + geom_errorbar(aes(xmin=boot_CI_lower, xmax=boot_CI_upper),position=position_dodge(width=0.7), width=.1,alpha=0.5,color='darkblue',show.legend = FALSE) 
    }## plot error bar or not
    print(plot)
  }
  return(output_df)}

plot_control_age_roc <- function(df, col=col_roc_E5, title=' ',age1="65", age2="75",age3="85", legendname=FALSE, APOE=FALSE, plot=TRUE){
  QC1 = control_age_roc(df,as.integer(age1),col=col)
  QC2 = control_age_roc(df,as.integer(age2),col=col)
  QC3 = control_age_roc(df,as.integer(age3),col=col)
  all = rbind(QC1, QC2, QC3) 
  all$age <- factor(all$age, levels = c(age1, age2, age3))
  all = process_prs_col_name(all)
  
  if(plot==TRUE){
    plot <- ggplot(data = all, aes(x=auc, y = PRS, shape = age))+
      geom_point(size=3,alpha=0.4,position = position_dodge(width = 0.7))+
      facet_wrap(~ethnicity, ncol=1)+ guides(col=guide_legend('control_age_thres'))+
      xlab('AUC')+ ggtitle(title)+xlim(min(min(all$boot_CI_lower),0.45),max(max(all$boot_CI_upper),0.75))+
      theme_bw() 
    plot <- plot + geom_errorbar(aes(xmin=boot_CI_lower, xmax=boot_CI_upper),position=position_dodge(width=0.7), width=.1,alpha=0.2,show.legend = FALSE) 
    print(plot)
    }
  return(all)
}

plot_control_age_roc_multi <- function(df1, df2, df3, df1_name, df2_name, df3_name, title, legendname=FALSE){
  df1$method =df1_name 
  df2$method =df2_name 
  df=rbind(df1,df2)
  if(class(df3) != "logical"){
    df3$method=df3_name
    df=rbind(df, df3)
    df$method = factor(df$method, level=c(df1_name, df2_name, df3_name))
  }else{
    df$method = factor(df$method, level=c(df1_name, df2_name))
  }
  plot <- ggplot(data = df, aes(x=auc, y = PRS, shape = age, color=method)) + scale_fill_hue()+
    geom_point(size=2,alpha=0.7,position = position_dodge(width = 0.7))+
    facet_wrap(~ethnicity, ncol=1)+
    xlab('AUC')+ xlim(0.42,0.72)+ggtitle(title)+guides(col=guide_legend("annotations"), shape=guide_legend('age threshold for controls'))+
    theme_bw() 
  plot <- plot + geom_errorbar(aes(xmin=boot_CI_lower, xmax=boot_CI_upper),position=position_dodge(width=0.7), width=.2,alpha=0.3,show.legend = FALSE) 
  if(legendname!=FALSE){
    plot=plot +
      guides(col=guide_legend(legendname))
  }
  return(plot)
}

plot_control_age_roc_multi(plot_control_age_roc(wightman_bl[-34], col=col_roc_polypred3, plot=FALSE),
                           plot_control_age_roc(wightman_enformer, col=col_roc_polypred3, plot=FALSE),
                           plot_control_age_roc(wightman_update_all, col=col_roc_polypred3, plot=FALSE),
                           'baseline','bl+enformer','all anno', title='Wightman')



wightman_polypred_control_age = plot_control_age_roc(wightman_polypred, col=col_roc_polypred3, title='All annotations')
wightman_susie_control_age = plot_control_age_roc(wightman_susie, col=col_roc_polypred3, title='no annotations', plot=FALSE)

plot_control_age_roc_multi(plot_control_age_roc(wightman_polypred_converged, col=col_roc_polypred3, plot=FALSE),
                           plot_control_age_roc(wightman_bl_converged, col=col_roc_polypred3, plot=FALSE),
                           plot_control_age_roc(wightman_susie_converged, col=col_roc_polypred3, plot=FALSE),
                           'all annotations', 'baseline','none', title='Wightman')

plot_control_age_roc_multi(plot_control_age_roc(kunkle_polypred, col=col_roc_polypred3, plot=FALSE),
                           plot_control_age_roc(kunkle_bl, col=col_roc_polypred3, plot=FALSE),
                           plot_control_age_roc(kunkle_susie, col=col_roc_polypred3, plot=FALSE),
                           'all','baseline','none','Kunkle')

plot_control_age_roc_multi(plot_control_age_roc(bellenguez_polypred, col=col_roc_polypred3, plot=FALSE),
                           plot_control_age_roc(bellenguez_bl, col=col_roc_polypred3, plot=FALSE),
                           plot_control_age_roc(bellenguez_susie, col=col_roc_polypred3, plot=FALSE),
                           'all','baseline','none','bellenguez')



plot_control_age_roc_multi(plot_control_age_roc(polypred_PRSCS, col=col_roc_E5, plot=FALSE),
                           plot_control_age_roc(polypred_PRSCS_NOT0, col=col_roc_E5, plot=FALSE),
                           plot_control_age_roc(polypred_prscs_noEnformer, col=col_roc_E5, plot=FALSE),
                          "polyfun_PRSCS(PIP >0.3)",  "polyfun_PRSCS(PIP > 0)", 'polyfun_PRSCS_noenformer(PIP>0)','genomewide_plink')

##others -----
sbayesr_log <-  glm(Diagnosis~PRS+Sex+Age,family = 'binomial', data = sbayesR)
prsice_log <- glm(Diagnosis~PRS+Sex+Age,family = 'binomial', data = PRSice)


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


## others
SNP4prscs <- read.csv('/gpfs/commons/home/tlin/output/wightman/new_anno_0203/update_all+enformer/finemap/max_snp_5/agg_extract_0.3.tsv',sep = '\t', header=T,fill = T)
SNP_prscsnot0 <- read.csv('/gpfs/commons/home/tlin/output/wightman/new_anno_0203/update_all+enformer/finemap/max_snp_5/agg_extract_pip_not0.tsv', sep = '\t',header=T,fill = T)

ggplot(SNP_prscsnot0, aes(x = P)) +
  geom_histogram(aes(y = ..density..),  color = "black", fill='lightblue',bins=50) +
  labs(title = "Wightman_SNP (PIP > 0)") +
  xlab("p value") + coord_flip()+  theme_bw()



ggplot(SNP_prscsnot0, aes(x = P)) + scale_x_log10()+
  geom_density(alpha=.2, fill="#FF6666")  # Overlay with transparent density plot

## 36k
lookup_table <- data.frame(Race = c(1, 2, 4, 5),
                           final_population = c("AMR", "ASN", "AFR", "EUR"))

# Merge data frames based on 'race'
#PRSCS_36K <- merge(PRSCS_36K, lookup_table, by = "Race", all.x = TRUE)
#plink_36k <- merge(plink_36k, lookup_table, by='Race', all.x=TRUE)
colnames(PRSCS_36k_new_ancestry) = gsub("_0\\.", "_", colnames(PRSCS_36k_new_ancestry))

plot_ethnic_roc(PRSCS_36k_new_ancestry, title='wightman_36k, new_ancestry', plot=TRUE)
plot_ethnic_roc(PRSCS_36K, title='wightman_36k', plot=TRUE)

plot_ethnic_roc(PRSCS_36k_fixed_hispanic, title='fix hispanic 36k PRSCS', plot=TRUE)

plot_ethnic_roc_facet(PRSCS_36K,PRSCS_36k_fixed_hispanic, PRSCS_36k_new_ancestry, legend = 'ancestry',title='wightman 36k PRSCS',  QC1name = 'original', QC2name ='fixed hispanic',QC3name ='1K genome')

ASN=roc_result_boot(extract_race(PRSCS_36K,'ASN'), column_for_roc = col_roc_E5, boot_num = 50)

y_ordered <- factor(y, levels = unique(y))
ASN=process_prs_col_name(ASN)

ASN <- ASN %>%
  mutate_at(vars(c('boot_mean','boot_CI_upper','boot_CI_lower')), as.numeric)

ggplot(data = ASN, aes(x=boot_mean, y = PRS))+
    geom_point(size=3,alpha=0.9,position = position_dodge(width = 0.7), color='darkblue')+
    xlab('AUC')+ ggtitle('ASN')+xlim(0.3, 0.72)+theme_bw() +
  geom_errorbar(aes(xmin=boot_CI_lower, xmax=boot_CI_upper),position=position_dodge(width=0.7), width=.1,alpha=0.5,color='darkblue',show.legend = FALSE) 

pie(table(PRSCS_36K$final_population))
ASN


PRSCS_36K_bellenguez <- merge(PRSCS_36K_bellenguez, lookup_table, by = "Race", all.x = TRUE)

plot_ethnic_roc(PRSCS_36K, title='bellenguez_36k', plot=TRUE, col=  list("PRS_01","PRS_05","PRS_1","PRS_5"))
ASN=roc_result_boot(extract_race(PRSCS_36K_bellenguez ,'ASN'), column_for_roc = list("PRS_01","PRS_05","PRS_1","PRS_5"), boot_num = 50)

ASN=process_prs_col_name(ASN)

ggplot(data = ASN, aes(x=boot_mean, y = PRS))+
  geom_point(size=3,alpha=0.9,position = position_dodge(width = 0.7), color='darkblue')+
  xlab('AUC')+ ggtitle('ASN')+xlim(0.3, 0.72)+theme_bw() +
  geom_errorbar(aes(xmin=boot_CI_lower, xmax=boot_CI_upper),position=position_dodge(width=0.7), width=.1,alpha=0.5,color='darkblue',show.legend = FALSE) 


