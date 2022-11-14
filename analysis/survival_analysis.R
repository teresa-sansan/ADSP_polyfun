install.packages(c("survival", "ggsurvfit", "gtsummary", "tidycmprsk"))
library(survival)
library(ggsurvfit)
#library(survminer)
library(gtsummary)
library(tidycmprsk)
#library(condsurv)
library('gridExtra')
library(grid)
#library(survminer)




par(mfrow=c(3,1))


surv_fit <- function(df, header){
   survfit2(Surv(Age, Diagnosis) ~ Sex , data = df) %>%
      ggsurvfit() + xlim(c(65,92)) +  add_confidence_interval(type='ribbon') +
      labs(x='age', y = 'Survival Probability') + ggtitle(header)+scale_color_manual(values = c("blue", "red"), labels = c("control", "case")) +
      scale_fill_manual(values = c("blue", "red"), labels = c("control", "case"))
}

cox_fit <- function(df){
   print(coxph(Surv(Age, Diagnosis) ~ Sex + PRS1 + X1 + X2 + X3 + X4 + X5+ X6+ X7+ X8+ X9+ X10, data = df))
}

plot_surv_ethnicity <- function(df, title_name){
   EUR <- surv_fit(extract_eur(df), 'EUR') + theme(legend.position="none")  + xlab(NULL)
   cox_fit(extract_eur(df))
   AFR <- surv_fit(extract_afr(df), 'AFR') + theme(legend.position="none")  + xlab(NULL)
   cox_fit(extract_afr(df))
   AMR <- surv_fit(extract_amr(df), 'AMR')
   cox_fit(extract_amr(df))
   title <- ggdraw() + draw_label(title_name, fontface='bold')
   plot <- plot_grid(EUR,AFR,AMR, nrow=3)
   plot_grid(title, plot, ncol=1, rel_heights=c(0.1, 1)) 
   
}


plot_surv_ethnicity(wightman_polypred, 'wightman_polypred (covariant: sex)')
plot_surv_ethnicity(kunkle_polypred, 'wightman_polypred (covariant: sex)')




