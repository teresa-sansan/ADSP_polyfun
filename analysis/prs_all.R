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
sbayesR = rename_prepreocess("/gpfs/commons/home/tlin/output/prs/sbayesR.tsv")
  
bellenguez_update <- rename_preprocess('/gpfs/commons/home/tlin/output/prs/bellenguez_bellenguez_updateRSID_prs_PC.tsv')
bellenguez_all2 <- pre_process("/gpfs/commons/home/tlin/output/prs/bellenguez_all_2/bellenguez_all_PC_CS.tsv")
bellenguez_cT <- pre_process('/gpfs/commons/home/tlin/output/prs/bellenguez/bellenguez_pT_PRS_withPC.tsv')

bellenguez_qc <- pre_process('/gpfs/commons/home/tlin/output/prs/bellenguez/bellenguez_qc_pT_PRS.tsv')
bellenguez_qc_on_target <- pre_process('/gpfs/commons/home/tlin/output/prs/bellenguez/qc_on_target.tsv')
bellenguez_qc_on_base <- pre_process('/gpfs/commons/home/tlin/output/prs/bellenguez/qc_on_base.tsv')

bellenguez_cT0224 <- pre_process('/gpfs/commons/home/tlin/output/prs/bellenguez/fixed_0224/bellenguez_no_qc.tsv')
bellenguez_qc0224 <- pre_process('/gpfs/commons/home/tlin/output/prs/bellenguez/fixed_0224/bellenguez_qc.tsv')
bellenguez_qc_on_target0224 <- pre_process('/gpfs/commons/home/tlin/output/prs/bellenguez/fixed_0224/bellenguez_qc_on_target.tsv')
bellenguez_qc_on_base0224 <- pre_process('/gpfs/commons/home/tlin/output/prs/bellenguez/fixed_0224/bellenguez_qc_on_base.tsv')

##susie
susie <- rename_preprocess("/gpfs/commons/home/tlin/output/prs/bellenguez_susie/bellenguez_susie_prs_PC.tsv")

##kunkle
kunkle_cT <- pre_process('/gpfs/commons/home/tlin/output/prs/kunkle/kunkle_pT_prs_before_qc.tsv')
only2 <- pre_process("/gpfs/commons/home/tlin/output/prs/kunkle/kunkle_2snp.tsv")
##kunkle_qc
kunkle_cT_qc <- pre_process('/gpfs/commons/home/tlin/output/prs/kunkle/kunkle_qc_pT_PRS.tsv')
kunkle_cT_qc_update <- pre_process('/gpfs/commons/home/tlin/output/prs/kunkle/kunkle_qc_pT_PRS_correct_cols.tsv')
kunkle_qc_target <-pre_process('/gpfs/commons/home/tlin/output/prs/kunkle/qc_on_target.tsv')
kunkle_qc_base <- pre_process('/gpfs/commons/home/tlin/output/prs/kunkle/qc_on_base.tsv')

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

roc_result <-function(df, title=' ',column_for_roc=col_roc){
  if('e5' %in% col_roc){
    roc_list <- roc(Diagnosis ~ PRS_e5+PRS_001+PRS_005+PRS_01+PRS_05+PRS_1+PRS_5, data = df)
  }else
    roc_list <- roc(Diagnosis ~ PRS_001+PRS_005+PRS_01+PRS_05+PRS_1+PRS_5, data = df)
  
  plot_roc(roc_list,column_for_roc, title=title)
}



###c+T (before and after qc)-----
#### bellenguez-----
##bellenguez
col_roc <- list("PRS_001","PRS_005","PRS_01","PRS_05","PRS_1","PRS_5")
bellenguez_cT_plot <- roc_result(bellenguez_cT,col_roc, title = "bellenguez, clumping+pT" )

##bellenguez_qc
bellenguez_qc_plot = roc_result(bellenguez_qc,col_roc, title = "bellenguez, c+pT, QC" )
plot_grid(bellenguez_cT_plot, bellenguez_qc_plot,ncol = 2, nrow = 1)

## bellenguez_EUR
bellenguez_eur_plot = roc_result(extract_eur(bellenguez_cT),col_roc, title = "bellenguez, c+pT, EUR" )

##bellenguez_qc_eur
bellenguez_qc_eur_plot = roc_result(extract_eur(bellenguez_qc),col_roc, title = "bellenguez, c+pT, QC, EUR" )

plot_grid(pt_plot, bellenguez_qc_plot,ncol = 2, nrow = 2)


####kunkle------

## Kunkle
### ori
kunkle_pt_plot = roc_result(kunkle_cT,"kunkle, clumping+pT")

### Qced
kunkle_pt_qc_plot_update=roc_result(kunkle_cT_qc_update,title="qc, clumping+pT")

### base
kunkle_pt_qc_plot_base=roc_result(kunkle_qc_base, title="kunkle_qc_summary_stats")

### target
kunkle_pt_qc_plot_target=roc_result(kunkle_qc_target,title="kunkle_qc_target")
plot_grid(kunkle_pt_plot,kunkle_pt_qc_plot_update,kunkle_pt_qc_plot_base,kunkle_pt_qc_plot_target,ncol = 2, nrow = 2)

## Kunkle_EUR
### ori
kunkle_pt_plot=roc_result( extract_eur(kunkle_cT),title="kunkle(EUR)")

### QCed EUR only
kunkle_pt_qc_plot_update=roc_result(extract_eur(kunkle_cT_qc_update),title="qc, clumping+pT (EUR)")

### base 
kunkle_pt_qc_plot_base=roc_result(extract_eur(kunkle_qc_base), title="kunkle_qc_summary_stats(EUR)")

### target
kunkle_pt_qc_plot_target=roc_result(extract_eur(kunkle_qc_target),title="kunkle_qc_target(EUR)")
plot_grid(kunkle_pt_plot,kunkle_pt_qc_plot_update,kunkle_pt_qc_plot_base,kunkle_pt_qc_plot_target,ncol = 2, nrow = 2)



### polyfun----

##bellenguez_all
roc_list <- roc(Diagnosis ~ PRS1+PRS3+PRS5+PRS10, data = bellenguez_all2)
col_roc <- list("PRS1","PRS3","PRS5",'PRS10')
bellenguez_roc = plot_roc(roc_list,col_roc, legend = 'PRS', replace = 'Max SNP = ', title="bellenguez(polyfun)")

##bellenguez_eur
eur_list <- roc(Diagnosis ~ PRS1+PRS3+PRS5+PRS10, data =bellenguez_eur)
eur_roc = plot_roc(eur_list, col_roc, legend = 'PRS', replace = 'Max SNP = ', title="bellenguez(EUR)")
plot_grid(bellenguez_roc, eur_roc,ncol = 2, nrow = 1)

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

log_reg <- function(df, prs,main_title,plot = TRUE, legend="PRS_", replace="pT = 0."){
  mod1 <- glm(Diagnosis ~ Sex + Age, data= df, family=binomial)
  mod1_R2 <- NagelkerkeR2(mod1)$R2
  for (i in 1:length(prs)){
    frm <- as.formula(paste("Diagnosis ~ Sex + Age+", prs[[i]]))
    mod <- glm(formula = frm,family = 'binomial', data = df)
    mod_R2 <- NagelkerkeR2(mod)$R2 
    ##plot
    if(plot == TRUE){
    plot(residuals(mod, type = "pearson") ~ df[,prs[[i]]], ylab = "Residual", xlab = "PRS", col = as.factor(df$Diagnosis), 
         main =  paste(str_replace(prs[[i]], legend,replace), ", in ", population, 'population'), cex = 0.5)
    plot(residuals(mod, type = "pearson") ~ residuals(mod1, type = "pearson"), ylab = "w. PRS", xlab = "w/o PRS", col = as.factor(df$Diagnosis),
         cex=0.5,main = paste(str_replace(prs[[i]], legend,replace), ", in ", population, 'population'))
    abline(0,1, col='blue')
    legend("bottomright", legend = c("control","case"), col = 1:2, pch = 19, bty = "n")
    mtext(paste("Diagnosis ~ Sex+Age+PRS, df=" ,main_title),                   # Add main title
          side = 3,
          line = -1.25,
          outer = TRUE, cex =1.1)
    }
    ## table
    if(i ==1){
      cof <- data.frame(round(summary(mod)$coefficients[,c(1,2,4)],4))
      cof['R2'] = round(mod1_R2,5)
      cof['R2_prs'] = round(mod_R2,5)
      cof['PRS'] =prs[i]
    }
    else{
      cof2 <- data.frame(round(summary(mod)$coefficients[,c(1,2,4)],4))
      cof2['R2'] = round(mod1_R2,5)
      cof2['R2_prs'] = round(mod_R2,5)
      cof2['PRS'] =prs[i]
      cof <- rbind(cof, cof2)
    }  
  }
  return(cof)
}

kunkle_cT_log<-log_reg(kunkle_cT, col_roc, "kunkle_cT","all", plot=FALSE)
bell_log <- log_reg(bellenguez_update, list("PRS1","PRS3","PRS5",'PRS10'), 'bellenguez_updated', plot = FALSE, legend="PRS", replace="max_snp_per_locus=")
bell_log_eur <- log_reg(bellenguez_eur, list("PRS1","PRS3","PRS5",'PRS10'), 'bellenguez_updated', plot=FALSE, legend="PRS", replace="max_snp_per_locus=")
bell_susie_log<- log_reg(susie, list("PRS1","PRS3","PRS5",'PRS10'), 'bellenguez_susie', 'all',plot = FALSE, legend="PRS", replace="max_snp_per_locus=")

bell_qc_log <- log_reg(bellenguez_cT_qc, col_roc, "bellenguez_cT (QCed)", plot=FALSE)
bell_qc_EUR_log <- log_reg(extract_eur(bellenguez_cT_qc), col_roc, "bellenguez_cT (QCed, EUR)", plot=FALSE)

kunkle_cT_log<-log_reg(extract_eur(kunkle_cT), col_roc, "kunkle_cT", plot=FALSE)
kunkle_cT_qc_log<-log_reg(extract_eur(kunkle_cT_qc), col_roc, "kunkle_QCed","all", plot=FALSE)
kunkle_cT_qc_target_log<-log_reg(extract_eur(kunkle_qc_target), col_roc, "kunkle_QCed","QC on target", plot=FALSE)
kunkle_cT_qc_base_log<-log_reg(extract_eur(kunkle_qc_base), col_roc, "kunkle_QCed","qx on summary stats", plot=FALSE)
kunkle_cT_qc_EUR_log<-log_reg(extract_eur(kunkle_cT_qc), col_roc, "kunkle_QCed (EUR)","all", plot=FALSE)

par(mfrow=c(2,2))
plotR2(kunkle_cT_log,"kunkle, cT, EUR")
plotR2(kunkle_cT_qc_log,"kunkle, qc, EUR")
plotR2(kunkle_cT_qc_target_log,"kunkle, qc on target, EUR")
plotR2(kunkle_cT_qc_base_log,"kunkle, qc on summary stats,EUR")


sbayesr_log <-  glm(Diagnosis~PRS+Sex+Age,family = 'binomial', data = sbayesR)
prsice_log <- glm(Diagnosis~PRS+Sex+Age,family = 'binomial', data = PRSice)

par(mfrow=c(1,1))
bell_cT <- log_reg(pT, c("PRS_005","PRS_05","PRS_5"), "bellenguez_cT","all")

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


plotR2 <- function(log_output, header){
  prs_col = seq(from=1, to=dim(log_output)[1], by=4)
  plot(log_output$R2_prs[prs_col], ylim  = c(min(log_output$R2[1], min(log_output$R2_prs[prs_col])), max(log_output$R2[1], max(log_output$R2_prs[prs_col]))), 
       ylab='R-squared', xlab = "different PRS",pch=4,col="darkgreen", main=header,xaxt = "n")
  text(seq(from=1,to=length(prs_col)),log_output$R2_prs[prs_col] , 
       labels = log_output$R2_prs[prs_col], adj = c(0.5,2), xpd = TRUE, cex=1) 
  
  
  abline(h=log_output$R2[1], col=c("blue"), lty=2, lwd=3)
  #text((length(prs_col)+1)/2,log_output$R2[1] , labels =  "Null model", adj = c(0.5,-2), xpd = TRUE, cex=1) 
  text((length(prs_col)+1)/2,log_output$R2[1] , labels =  paste("Null model: ",log_output$R2[1]), adj = c(0.5,-1.5), xpd = TRUE, cex=1)
  axis(1,                         # Define x-axis manually
       at = 1:length(prs_col),
       labels = c(str_replace(log_output$PRS[1],"PRS_e5","p = 1*e-5"),str_replace(log_output$PRS[seq(from=1, to=dim(test)[1], by=4)][-1],"PRS_","p = 0.")))
}

plotR2(kunkle_cT_qc_EUR_log,"kunkle, qc, EUR")
plotR2(kunkle_cT_qc_log,"kunkle, qc")


