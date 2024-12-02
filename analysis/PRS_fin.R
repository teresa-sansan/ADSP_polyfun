library(stringr)
library(tidyr)
library(patchwork)
library(pROC)
library(readr)
library(dplyr)

## plot snp count ------------
# snp_count = read.table("/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/ADSP_reference_panel/fine_mapping/annotations_dl/finemap_v3_backup_teresa_update/pip_thres/credible_set_snp.count")
# colnames(snp_count) <- c('anno','pip_thres','snp_counts_one', 'snp_counts_all')
# snp_count$added_snp = snp_count$snp_counts_all- snp_count$snp_counts_one
#snp_count <- snp_count[snp_count$pip_thres != 0.0,]

snp_count = read.table('/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/PRS/36k_hg38/snp_pip_thres/bellenguez_pip_count.txt')
colnames(snp_count) <- c('anno','pip_thres','snp_counts_one')
snp_count$anno <- factor(snp_count$anno, levels = c('susie', 'baseline', 'omics', 'omics_dl'))

snp_count %>%
  ggplot(aes(x = pip_thres , y = snp_counts_one, fill = anno)) +
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  geom_text(aes(label = snp_counts_one, vjust = -0.3), position = position_dodge(width = 0.09), size = 3)+  
  xlab('pip threshold ')  + ylab('snp') + 
  scale_x_continuous(breaks = c(seq(from = 0, to = 0.9, by = 0.1)))+
  theme_minimal()+
  ggtitle('SNP counts with PIP threshold')

plot1
# plot2
# plot2<-snp_count %>%
#   ggplot(aes(x = pip_thres , y = snp_counts_all, fill = anno)) +
#   geom_bar(stat="identity", color="black", position=position_dodge()) +
#   geom_text(aes(label = snp_counts_all, vjust = -0.3), position = position_dodge(width = 0.09), size = 2.5)+  
#   xlab('pip threshold')  + ylab('snp') + 
#   scale_x_continuous(breaks = c(seq(from = 0, to = 0.9, by = 0.1)))+
#   theme_minimal()+
#   ggtitle('SNP counts (credible set based)')
# 
# plot1/plot2

## update snp count



## load prs (remove index0, aug)------------
path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/ADSP_reference_panel/fine_mapping/annotations_dl/finemap_v3_backup_teresa_update/prs/'
bellenguez_adsp_susie_aug <- pre_process(paste(path, 'prs_susie.tsv', sep = ''))
bellenguez_adsp_baseline_aug <- pre_process(paste(path, 'prs_baseline.tsv', sep = ''))
bellenguez_adsp_omics_aug <- pre_process(paste(path, 'prs_omics.tsv', sep = ''))
bellenguez_adsp_omics_dl_aug <- pre_process(paste(path, 'prs_omics_dl.tsv', sep = ''))

## load prs (remove index0, fin)------
path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/PRS/36k_hg38/polyfun/pip_thres/'
bellenguez_adsp_susie <- pre_process(paste(path, 'susie.tsv', sep = ''))
bellenguez_adsp_baseline <- pre_process(paste(path, 'baseline.tsv', sep = ''))
bellenguez_adsp_omics <- pre_process(paste(path, 'omics.tsv', sep = ''))
bellenguez_adsp_omics_dl <- pre_process(paste(path, 'omics_dl.tsv', sep = ''))

## c+ PT
Cpt <- pre_process('/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/PRS/36k_hg38/c+pt/c_pt.tsv')
p_list = list("P_e5","P_e4","P_0.001","P_0.01","P_0.1", "P_0.0" )
plot_ethnic_roc(Cpt, col =p_list, title='c+pt', plot=TRUE)

##prscs
#bellenguez_prscs <- pre_process('/gpfs/commons/home/tlin/output/prs/PRSCS/36k_adsp_ld_panel/bellenguez/_bellenguez_adsp_ld_36k.tsv')
prscs <- pre_process('/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/PRS/36k_hg38/prscs/bellenguez_adsp_ld_36k.tsv')
prscs_ukbb <- pre_process('/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/PRS/36k_hg38/prscs/ukbb/bellenguez_ukbb_ld_36k.tsv')
plot_ethnic_roc(prscs_ukbb, col =list("P_e5","P_e4","P_0.001","P_0.01","P_0.1" ,'P_0.0'), title='prscs_ukbb', plot=TRUE)
plot_ethnic_roc(prscs, col =list("P_e5","P_e4","P_0.001","P_0.01","P_0.1" ,'P_0.0'), title='prscs', plot=TRUE)


## ROC result
cpt_prscs_roc = plot_ethnic_roc_facet( Cpt,prscs, prscs_ukbb, QC2name="PRSCS_ADSP", QC3name="PRSCS_UKBB",QC1name="C+pT", 
                                      col = list("P_e5","P_e4","P_0.001","P_0.01","P_0.1" ,'P_0.0'), title = '' ,  legendname = 'method', return_plot=FALSE)
plot_ethnic_roc_facet(prscs, prscs_ukbb,  data.frame(),QC1name="PRSCS_ADSP", QC2name="PRSCS_UKBB", 
                       col = list("P_e5","P_e4","P_0.001","P_0.01","P_0.1" ,'P_0.0'), title = '' ,  legendname = 'method', return_plot=TRUE)
# cpt_prscs_roc = plot_ethnic_roc_facet(Cpt, prscs, data.frame(), QC1name="C+pT", QC2name="PRSCS_ADSP",QC3name="omics", 
#                       col = list("P_e5","P_e4","P_0.001","P_0.01","P_0.1" ,'P_0.0'), title = '' ,  legendname = 'method')

col_roc_pthres <- list("PRS_0","PRS_1","PRS_2","PRS_3","PRS_4","PRS_5","PRS_6","PRS_7","PRS_8", "PRS_9")
polyfun_roc = plot_ethnic_roc_facet(bellenguez_adsp_susie, bellenguez_adsp_baseline, bellenguez_adsp_omics, bellenguez_adsp_omics_dl, 
                                    QC1name="SuSiE", QC2name="bl",QC3name="bl/ omics", QC4name="bl/ omics/ dl",  
                                    col = col_roc_pthres, title = 'remove index 0' , legendname = 'annotation', return_plot = FALSE)

roc = rbind(polyfun_roc, cpt_prscs_roc)
ggplot(roc, aes(x = qc_status, y = auc, fill = qc_status)) +
  geom_boxplot(alpha = 0.8, outlier.size = 1) +
  facet_wrap(~ethnicity, ncol = 3) +
  labs(
    x = "",
    y = "auROC",
    fill = "PRS"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1),
    panel.grid.major = element_line(color = "gray90"),
    strip.background = element_rect(fill = "gray90"),
    strip.text = element_text(face = "bold"),legend.position = "none"
  )



col_roc_pthres_prscs <- list("P_0.0001","P_0.001","P_0.1")
roc_result(bellenguez_prscs, column_for_roc = col_roc_pthres_prscs)

roc_result <-function(df, column_for_roc = col_roc_E5){
  auc_list=list()
  auc_df = data.frame("PRS" = unlist(column_for_roc), "auc" = 0)
  print(auc_df)
  for (i in 1:length(column_for_roc)){
    col = column_for_roc[[i]]
    print(col)
    #roc_formula = roc(df[["Diagnosis"]],df[[col]], quiet=T)
    print(df[["Diagnosis"]])
    print(df[[col]])
    auc_value <- roc(df[["Diagnosis"]],df[[col]] , quiet=T)$auc
    auc_list[i] = as.numeric(auc_value)
  }
  auc_df$auc = unlist(auc_list)
  return(auc_df)
}

## get beta for logistic regression
get_beta <- function(df = bellenguez_adsp_susie, col = col_roc_pthres){
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

relative_r2 <- function(df, prs) {
  # Helper function for Nagelkerke's R²
  nagelkerke_R2 <- function(LL_model, LL_null, n) {
    1 - (exp(LL_null - LL_model)^(2 / n)) / (1 - exp(2 * LL_null / n))
  }
  
  ancestries <- c("EUR", "AFR", "AMR")
  results <- data.frame(
    Ancestry = character(),
    PRS = character(),
    R2_base = numeric(),
    R2_full = numeric(),
    Relative_R2 = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Loop through ancestries
  for (ancestry in ancestries) {
    df_ancestry <- df[df$predicted_ancestry == ancestry, ]
    
    # Fit the base model (covariates only)
    mod_base <- glm(Diagnosis ~ Sex + Age, data = df_ancestry, family = binomial)
    LL_base <- logLik(mod_base)
    LL_null <- logLik(glm(Diagnosis ~ 1, data = df_ancestry, family = binomial))  # Null model
    n <- nrow(df_ancestry)
    
    # Loop through each PRS variable
    for (i in seq_along(prs)) {
      # Fit the full model (covariates + PRS)
      formula_full <- as.formula(paste("Diagnosis ~ Sex + Age +", prs[i]))
      mod_full <- glm(formula_full, data = df_ancestry, family = binomial)
      LL_full <- logLik(mod_full)
      
      # Calculate Nagelkerke's R²
      R2_base <- nagelkerke_R2(LL_base, LL_null, n)
      R2_full <- nagelkerke_R2(LL_full, LL_null, n)
      relative_R2 <- R2_full - R2_base
      
      # Store results
      results <- rbind(
        results,
        data.frame(
          Ancestry = ancestry,
          PRS = prs[i],
          R2_base = R2_base,
          R2_full = R2_full,
          Relative_R2 = relative_R2
        )
      )
    }
  }
  
  return(results)
}

df_r2 <- function(df, col, df_name){
  df_r2=relative_r2(df, col)
  df_r2$data = df_name
  return(df_r2)
}

r2=rbind(df_r2(prscs, c("P_e5","P_e4","P_0.001","P_0.01","P_0.1" ), 'PRSCS_adsp'),
         df_r2(prscs_ukbb, c("P_e5","P_e4","P_0.001","P_0.01","P_0.1" ), 'PRSCS_ukbb'),
         df_r2(Cpt, c("P_e5","P_e4","P_0.001","P_0.01","P_0.1" ), 'c_pT'), 
         df_r2(bellenguez_adsp_susie, unlist(col_roc_pthres), 'SuSiE'),
         df_r2(bellenguez_adsp_susie, unlist(col_roc_pthres), 'bl'),
         df_r2(bellenguez_adsp_susie, unlist(col_roc_pthres), 'bl/ omics'),
         df_r2(bellenguez_adsp_susie, unlist(col_roc_pthres), 'bl/ omics/ dl'))

r2$Ancestry <- factor(r2$Ancestry, levels = c('EUR', 'AFR', 'AMR'))
r2$data <- factor(r2$data, levels = c('SuSiE', 'bl', 'bl/ omics', 'bl/ omics/ dl', 'c_pT','PRSCS_adsp', 'PRSCS_ukbb'))


library(cowplot)
roc_plot <- ggplot(roc, aes(x = qc_status, y = auc, fill = qc_status)) +
  geom_boxplot(alpha = 0.8, outlier.size = 1) +
  facet_wrap(~ethnicity, ncol = 3) +
  labs(
    x = "",
    y = "auROC",
    fill = "PRS"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1),
    panel.grid.major = element_line(color = "gray90"),
    strip.background = element_rect(fill = "gray90"),
    strip.text = element_text(face = "bold"),legend.position = "none"
  )



r2_plot <-ggplot(data = r2, aes(x = data, y = Relative_R2, fill = data)) +
  geom_boxplot(alpha = 0.8, outlier.size = 1) +
  facet_wrap(~Ancestry, ncol = 3) +
  labs(
    x = "",
    y = "R²",
    fill = "PRS"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1),
    panel.grid.major = element_line(color = "gray90"),
    strip.background = element_rect(fill = "gray90"),
    strip.text = element_text(face = "bold"),legend.position = "none"
  )


combined_plot <- plot_grid(
  roc_plot,
  r2_plot,
  ncol = 1            
)
print(combined_plot)

