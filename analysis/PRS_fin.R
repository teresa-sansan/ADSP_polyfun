library(stringr)
library(tidyr)
library(patchwork)
library(pROC)

library(dplyr)

## plot snp count ------------
snp_count = read.table("/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/ADSP_reference_panel/fine_mapping/annotations_dl/finemap_v3_backup_teresa_update/pip_thres/credible_set_snp.count")
colnames(snp_count) <- c('anno','pip_thres','snp_counts_one', 'snp_counts_all')
snp_count$added_snp = snp_count$snp_counts_all- snp_count$snp_counts_one

snp_count <- snp_count[snp_count$pip_thres != 0.0,]
plot1 <- snp_count %>%
  ggplot(aes(x = pip_thres , y = snp_counts_one, fill = anno)) +
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  geom_text(aes(label = snp_counts_one, vjust = -0.3), position = position_dodge(width = 0.09), size = 2.5)+  
  xlab('pip threshold ')  + ylab('snp') + 
  scale_x_continuous(breaks = c(seq(from = 0, to = 0.9, by = 0.1)))+
  theme_minimal()+
  ggtitle('SNP counts (snp based)')


plot2<-snp_count %>%
  ggplot(aes(x = pip_thres , y = snp_counts_all, fill = anno)) +
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  geom_text(aes(label = snp_counts_all, vjust = -0.3), position = position_dodge(width = 0.09), size = 2.5)+  
  xlab('pip threshold')  + ylab('snp') + 
  scale_x_continuous(breaks = c(seq(from = 0, to = 0.9, by = 0.1)))+
  theme_minimal()+
  ggtitle('SNP counts (credible set based)')

plot1/plot2


## load prs (remove index0)------------
path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/ADSP_reference_panel/fine_mapping/annotations_dl/finemap_v3_backup_teresa/prs/'
bellenguez_adsp_susie_aug <- pre_process(paste(path, 'prs_susie.tsv', sep = ''))
bellenguez_adsp_baseline_aug <- pre_process(paste(path, 'prs_baseline.tsv', sep = ''))
bellenguez_adsp_omics_aug <- pre_process(paste(path, 'prs_omics.tsv', sep = ''))
bellenguez_adsp_omics_dl_aug <- pre_process(paste(path, 'prs_omics_dl.tsv', sep = ''))

bellenguez_adsp_susie_aug_noqc <- pre_process(paste(path, 'prs_no_qc_susie.tsv', sep = ''))
bellenguez_adsp_baseline_aug_noqc <- pre_process(paste(path, 'prs_no_qc_baseline.tsv', sep = ''))
bellenguez_adsp_omics_aug_noqc <- pre_process(paste(path, 'prs_no_qc_omics.tsv', sep = ''))
bellenguez_adsp_omics_dl_aug_noqc <- pre_process(paste(path, 'prs_no_qc_omics_dl.tsv', sep = ''))

col_roc_pthres <- list("PRS_0","PRS_1","PRS_2","PRS_3","PRS_4","PRS_5","PRS_6","PRS_7","PRS_8", "PRS_9")

check <-plot_ethnic_roc_facet(bellenguez_adsp_susie_aug, bellenguez_adsp_baseline_aug, bellenguez_adsp_omics_aug, bellenguez_adsp_omics_dl_aug, QC1name="susie", QC2name="baseline",QC3name="omics", 
                      QC4name="omics_dl",  col = col_roc_pthres, title = 'remove index 0' , boot_num=1, legendname = 'annotation')

#plot_ethnic_roc_facet(bellenguez_adsp_susie_aug_noqc, bellenguez_adsp_baseline_aug_noqc,bellenguez_adsp_omics_aug_noqc, bellenguez_adsp_omics_dl_aug_noqc, QC1name="susie", QC2name="baseline",QC3name="omics", 
#                      QC4name="omics_dl",  col = col_roc_pthres, title = 'remove index 0, no qc on genotype data' , boot_num=5, legendname = 'annotation')

pdf("/gpfs/commons/home/tlin/pic/prs.pdf",    width = 10,    height = 8)
plot_ethnic_roc_facet(bellenguez_adsp_susie_aug, bellenguez_adsp_baseline_aug, bellenguez_adsp_omics_aug, bellenguez_adsp_omics_dl_aug, QC1name="susie", QC2name="baseline",QC3name="omics", 
                       QC4name="omics_dl",  col = col_roc_pthres, title = 'remove index 0' , boot_num=50, legendname = 'annotation')

dev.off() 




roc_result <-function(df, column_for_roc = col_roc_E5){
  auc_list=list()
  auc_df = data.frame("PRS" = unlist(column_for_roc), "auc" = 0)
  for (i in 1:length(column_for_roc)){
    col = column_for_roc[[i]]
    #roc_formula = roc(df[["Diagnosis"]],df[[col]], quiet=T)
    auc_value <- roc(df[["Diagnosis"]],df[[col]] , quiet=T)$auc
    
    auc_list[i] = as.numeric(auc_value)
  }
  auc_df$auc = unlist(auc_list)
  return(auc_df)
}


get_beta <- function(df = bellenguez_adsp_susie_aug, col = col_roc_pthres){
   df<- df %>% mutate_at(unlist(col), ~(scale(.) %>% as.vector))
   beta <- list()
   beta_sd <- list()
   for (x in col) {
     model <- glm(as.formula(paste("Diagnosis ~", x)), data = df, family = binomial)
     beta[[x]] <- summary(model)$coefficients[x, "Estimate"]
     beta_sd[[x]] <- summary(model)$coefficients[x, "Std. Error"]
   }
   beta <- unlist(beta)
   beta_sd <- unlist(beta_sd)
   print(beta)
   print(beta_sd)
}

get_beta()
model <- glm(Diagnosis ~ PRS_0, data = bellenguez_adsp_susie_aug, family = binomial)
summary_model <- summary(model)
confint(model)

# Extract the beta coefficient for PRS_0
beta_PRS_0 <- summary_model$coefficients["PRS_0", "Estimate"]

# Print the beta coefficient
print(beta_PRS_0)



get_beta()
  