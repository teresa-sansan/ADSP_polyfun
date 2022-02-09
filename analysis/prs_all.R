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
sbayesR = rename_prepreocess("/gpfs/commons/home/tlin/output/prs/sbayesR.tsv")
  
bellenguez_update <- rename_preprocess('/gpfs/commons/home/tlin/output/prs/bellenguez_updateRSID_prs_PC.tsv')
bellenguez_all2 <- pre_process("/gpfs/commons/home/tlin/output/prs/bellenguez_all_2/bellenguez_all_PC_CS.tsv")

##susie
susie <- rename_preprocess(read.delim("/gpfs/commons/home/tlin/output/prs/bellenguez_susie/bellenguez_susie_prs_PC.tsv"))

##kunkle
kunkle_cT <- rename_preprocess('/gpfs/commons/home/tlin/output/cT/kunkle/kunkle_prs.tsv')


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

##AUC ------
plot_roc <- function(roclist, roccol,legend="PRS_", replace="p = 0.",title = FALSE){
  legendname = list()
  for (i in roccol){
    legendname = append(legendname,paste(str_replace(i, legend,replace),'auc =',round(roc_list[[i]]$auc,3)))
  }
  print(legendname)
  curve <- ggroc(roclist)
  curve + xlab("FPR") + ylab("TPR") + ggtitle(title)+theme_bw()+
    scale_colour_discrete(name="PRS",labels = c(legendname))
}   

##c+T
roc_list <- roc(Diagnosis ~ PRS_005+PRS_05+PRS_5, data = pT)
col_roc <- list("PRS_005","PRS_05","PRS_5")
pt_plot = plot_roc(roc_list,col_roc, title="bellenguez, clumping+pT")

roc_list <- roc(Diagnosis ~ PRS_005+PRS_05+PRS_5, data = kunkle_cT)
kunkle_pt_plot = plot_roc(roc_list,col_roc, title="kunkle, clumping+pT")

plot_grid(pt_plot, kunkle_pt_plot,ncol = 2, nrow = 1)

##bellenguez_all
roc_list <- roc(Diagnosis ~ PRS1+PRS3+PRS5+PRS10, data = bellenguez_all2)
col_roc <- list("PRS1","PRS3","PRS5",'PRS10')
bellenguez_roc = plot_roc(roc_list,col_roc, legend = 'PRS', replace = 'Max SNP = ', title="bellenguez(polyfun)")

##susie
roc_list <- roc(Diagnosis ~ PRS1+PRS3+PRS5+PRS10, data = susie)
susie_roc = plot_roc(roc_list,col_roc, legend = 'PRS', replace = 'Max SNP = ', title="bellenguez(SuSiE)")
plot_grid(bellenguez_roc,susie_roc,ncol = 2, nrow = 1)


## Regression ----

log_reg <- function(df, prs,main_title, population, population_subset = FALSE,plot = TRUE, legend="PRS_", replace="pT = 0."){
  if(population_subset == TRUE){
    df = df[df$final_population == population,]
  }
  mod1 <- glm(Diagnosis ~ Sex + Age, data= df, family=binomial)
  for (i in 1:length(prs)){
    frm <- as.formula(paste("Diagnosis ~ Sex + Age+", prs[[i]]))
    mod <- glm(formula = frm,family = 'binomial', data = df)
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
      cof['PRS'] =prs[i]
    }
    else{
      cof2 <- data.frame(round(summary(mod)$coefficients[,c(1,2,4)],4))
      cof2['PRS'] =prs[i]
      cof <- rbind(cof, cof2)
    }  
  }
  return(cof)
}
log_reg(kunkle_cT, col_roc, "kunkle_cT","all")
bell_log <- log_reg(bellenguez_update, list("PRS1","PRS3","PRS5",'PRS10'), 'bellenguez_updated', 'all',plot = FALSE, legend="PRS", replace="max_snp_per_locus=")
bell_log_eur <- log_reg(bellenguez_update, list("PRS1","PRS3","PRS5",'PRS10'), 'bellenguez_updated', 'EUR',population_subset = TRUE, legend="PRS", replace="max_snp_per_locus=")
par(mfrow=c(1,2))
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
p_beta_plot(bell_log,c(4,8,12,16),'bellenguez_updated')
p_beta_plot(bell_log_eur,c(4,8,12,16),'bellenguez_updated (EUR only)')

