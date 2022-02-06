library(ggplot2)
library(lattice)
library(dplyr)
library(stringr)
library(tidyr)
library(ggstance)
library(jtools)
library(tibble)


## A function that filtered out all the young samples and the samples without Diagnosis
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






# LOAD DATA ---------------
polypred1 = read.table("/gpfs/commons/home/tlin/output/prs/bellenguez_all_2/prs_diagnosis_0219.2021_max_snp_1.tsv",header=T,fill = T)
polypred3 = read.table("/gpfs/commons/home/tlin/output/prs/bellenguez_all_2/prs_diagnosis_0219.2021_max_snp_3.tsv",header=T,fill = T)
polypred5 = read.table("/gpfs/commons/home/tlin/output/prs/bellenguez_all_2/prs_diagnosis_0219.2021_max_snp_5.tsv",header=T,fill = T)
polypred7 = read.table("/gpfs/commons/home/tlin/output/prs/bellenguez_all_2/prs_diagnosis_0219.2021_max_snp_7.tsv",header=T,fill = T)
polypred10 = read.table("/gpfs/commons/home/tlin/output/prs/bellenguez_all_2/prs_diagnosis_0219.2021_max_snp_10.tsv",header=T,fill = T)


susie_polypred1 = read.table("/gpfs/commons/home/tlin/output/prs/bellenguez_susie/susie_prs_diagnosis_0219.2021_max_snp_1.tsv",header=T,fill = T)
susie_polypred3 = read.table("/gpfs/commons/home/tlin/output/prs/bellenguez_susie/susie_prs_diagnosis_0219.2021_max_snp_3.tsv",header=T,fill = T)
susie_polypred5 = read.table("/gpfs/commons/home/tlin/output/prs/bellenguez_susie/susie_prs_diagnosis_0219.2021_max_snp_5.tsv",header=T,fill = T)
susie_polypred7 = read.table("/gpfs/commons/home/tlin/output/prs/bellenguez_susie/susie_prs_diagnosis_0219.2021_max_snp_7.tsv",header=T,fill = T)
susie_polypred10 = read.table("/gpfs/commons/home/tlin/output/prs/bellenguez_susie/susie_prs_diagnosis_0219.2021_max_snp_10.tsv",header=T,fill = T)


beta_polypred1 = read.table("/gpfs/commons/home/tlin/output/prs/bellenguez_all_2/adjbeta_prs_diagnosis_0219.2021_max_snp_1.tsv",header=T,fill = T)
beta_polypred3 = read.table("/gpfs/commons/home/tlin/output/prs/bellenguez_all_2/adjbeta_prs_diagnosis_0219.2021_max_snp_3.tsv",header=T,fill = T)
beta_polypred5 = read.table("/gpfs/commons/home/tlin/output/prs/bellenguez_all_2/adjbeta_prs_diagnosis_0219.2021_max_snp_5.tsv",header=T,fill = T)
beta_polypred7 = read.table("/gpfs/commons/home/tlin/output/prs/bellenguez_all_2/adjbeta_prs_diagnosis_0219.2021_max_snp_7.tsv",header=T,fill = T)
beta_polypred10 = read.table("/gpfs/commons/home/tlin/output/prs/bellenguez_all_2/adjbeta_prs_diagnosis_0219.2021_max_snp_10.tsv",header=T,fill = T)

updatePRS1 = pre_process("/gpfs/commons/home/tlin/output/prs/bellenguez_all_2/with_PC/UPDATEprs_diagnosis_0219.2021_max_snp_1_subset.tsv")
updatePRS3 = pre_process("/gpfs/commons/home/tlin/output/prs/bellenguez_all_2/with_PC/UPDATEprs_diagnosis_0219.2021_max_snp_3_subset.tsv")
updatePRS5 = pre_process("/gpfs/commons/home/tlin/output/prs/bellenguez_all_2/with_PC/UPDATEprs_diagnosis_0219.2021_max_snp_5_subset.tsv")
updatePRS7 = pre_process("/gpfs/commons/home/tlin/output/prs/bellenguez_all_2/with_PC/UPDATEprs_diagnosis_0219.2021_max_snp_7_subset.tsv")
updatePRS10 = pre_process("/gpfs/commons/home/tlin/output/prs/bellenguez_all_2/with_PC/UPDATEprs_diagnosis_0219.2021_max_snp_10_subset_new.tsv")



## bind df and filter
bind_filt <- function(df1, df2, df3, df4, df5, preprocess = FALSE){
  df <- rbind(df1, df2, df3, df4, df5)
  snp = c(1,3,5,7,10)
  if(dim(df1)[1] == dim(df2)[1])
    df$SNP = rep(snp,each=dim(df1)[1])
  if(preprocess ==TRUE){
    df = pre_process(df, FILE = TRUE)
  return(df)
  }
  
}
polypred_filt = bind_filt(polypred1, polypred3, polypred5, polypred7, polypred10, preprocess = TRUE)
susie_polypred = rbind(susie_polypred1, susie_polypred3, susie_polypred5, susie_polypred7, susie_polypred10)
beta_polypred = rbind(beta_polypred1, beta_polypred3, beta_polypred5,beta_polypred7, beta_polypred10)

# --------
# polypred = rbind(polypred1, polypred3, polypred5, polypred7, polypred10)
# susie_polypred = rbind(susie_polypred1, susie_polypred3, susie_polypred5, susie_polypred7, susie_polypred10)
# beta_polypred = rbind(beta_polypred1, beta_polypred3, beta_polypred5,beta_polypred7, beta_polypred10)
# 
# snp = c(1,3,5,7,10)
# polypred$SNP = rep(snp,each=dim(polypred1)[1])
# susie_polypred$SNP = rep(snp,each=dim(susie_polypred1)[1])
# beta_polypred$SNP = rep(snp,each=dim(beta_polypred1)[1])
# 
# polypred <- subset(polypred, polypred$Diagnosis!= -1) 
# susie_polypred <- subset(susie_polypred, susie_polypred$Diagnosis!= -1) 
# beta_polypred <- subset(beta_polypred, beta_polypred$Diagnosis!= -1)


polypred$Diagnosis_case = "Control"
polypred[polypred$Diagnosis == 1,]$diagnosis_case = "Case"

susie_polypred$Diagnosis_case = "Control"
susie_polypred[susie_polypred$Diagnosis == 1,]$Diagnosis_case = "Case"


beta_polypred$Diagnosis_case = "Control"
beta_polypred[beta_polypred$Diagnosis == 1,]$Diagnosis_case = "Case"

sort_prs <-polypred[order(PRS),]
head(sort_prs)

draw_plot <- function(PRS, percentage,main, data){
  par(mfrow=c(2,1))
  num <- round(dim(sort_prs)[1]*percentage*0.01)
  plot(head(sort_prs, num)$PRS, head(sort_prs,num)$Diagnosis,ylim = c(0,1), xlab="PRS",ylab="Diagnosis",sub = paste("top ",percentage, "%"), main=main , data = data)
  plot(tail(sort_prs, num)$PRS, tail(sort_prs,num)$Diagnosis,ylim =c(0,1),xlab="PRS",ylab="diagnosis",sub = paste("bottom ",percentage, "%"))
}

draw_plot(polypred1,10, main = "SNP = 1", data = polypred)
draw_plot(polypred1,10, main = "SNP = 10")
draw_plot(polypred, 5) ## all

control = polypred[polypred$diagnosis == 0,]$PRS
case = polypred[polypred$diagnosis == 1,]$PRS


par(mfrow=c(2,1))
hist(polypred[polypred$Diagnosis == 1,]$PRS, main = "AD patients",xlab="PRS", ylab='number', freq= FALSE)
hist(polypred[polypred$Diagnosis == 0,]$PRS, main = "AD controls",xlab="PRS", ylab='number', freq= FALSE)

par(mfrow=c(2,1))
densityplot(~ polypred$PRS, polypred[polypred$diagnosis == 1,], main = "Case")
densityplot(~ polypred$PRS, polypred[polypred$diagnosis == 0,], main = "Control")

#densityplot(~PRS, data = polypred, group = diagnosis_case, auto.key= TRUE, plot.points=FALSE)


ggplot(polypred[polypred$diagnosis == 0,]) +
  geom_histogram(aes(x = PRS,y=..density..),
                 fill = "grey", color="black", bins=30)
  #geom_line(data =  aes(x = x, y = y), color = "red"))


snp.f <- factor(polypred$SNP, levels= c(1,3,5,7,10),
                labels = c("1 snp", "3 snp", "5 snp", "7 snp", "10 snp" ))

sm.density.compare(polypred$PRS, snp.f)
colfill<-c(2:(2+length(levels(snp.f))))
legend('topright', levels(snp.f),fill=colfill)


diagnosis.f <- factor(polypred$Diagnosis, levels= c(1,0),
                labels = c("case", "control"))

sm.density.compare(polypred1$PRS, diagnosis.f)
colfill<-c(2:(2+length(levels(diagnosis.f))))
legend('topright', levels(diagnosis.f),fill=colfill)


density_plot <- function(polypred, name){
  sort = polypred[order(polypred$PRS),]$PRS
  plot(density(polypred[polypred$Diagnosis == 1,]$PRS),col = "red", main=name, xlab="PRS") 
  lines(density(polypred[polypred$Diagnosis == 0,]$PRS), col = "blue")
  legend("topright", legend=c("Case", "Control"),
         col=c("red", "blue"), lty=1:1, cex=0.8,
         box.lty=0)
}

density_plot(polypred1, "snp per locus = 1")
legend("topleft", legend=c("Case", "Control"),
       col=c("red", "blue"), lty=1:1, cex=0.8,
       box.lty=0)

par(mfrow=c(2,1))
density_plot(polypred10, "snp per locus = 10 (OLD)")
legend("topleft", legend=c("Case", "Control"),
       col=c("red", "blue"), lty=1:1, cex=0.8,
       box.lty=0)

density_plot(updatePRS10, "snp per locus = 10 (UPDATE)")
legend("topleft", legend=c("Case", "Control"),
       col=c("red", "blue"), lty=1:1, cex=0.8,
       box.lty=0)


#TEST OUT ----

par(mfrow=c(2,1))
draw_top_densityplot <- function(polypred, percentage, main){
  sort_prs <-polypred[order(polypred$PRS),]
  num <- round(dim(sort_prs)[1]*percentage*0.01)
  head <- head(sort_prs,num)
  tail <- tail(sort_prs,num)
  plot(density(head$PRS, head$Diagnosis), xlab="PRS",sub = paste("top ",percentage, "%"), main = main)
  plot(density(tail$PRS, tail$Diagnosis), xlab="PRS",sub = paste("bottom ",percentage, "%"), main =""  )
}

draw_top_densityplot(polypred1, 10, "SNP1")


par(mfrow=c(2,2))



density_plot(polypred3, "snp per locus = 3")
density_plot(polypred5, "snp per locus = 5")
density_plot(polypred7, "snp per locus = 7")
density_plot(polypred10, "snp per locus = 10")

legend("topleft", legend=c("Case", "Control"),
       col=c("red", "blue"), lty=1:1, cex=0.8,
       box.lty=0) 



# Subseting >2SD
polypred_sd = polypred[abs(scale(polypred$PRS))>2,]
susie_polypred_sd = susie_polypred[abs(scale(susie_polypred$PRS))>2,]


#top and tail %
percentage_subset<- function (df, percent){
  num = round(dim(df)[1]*percent*0.01)
  df = df[order(df$PRS),]
  data1 = head(df, num)
  data2 = tail(df, num)
  return(rbind(data1, data2))
}



# check out how simlar is SuSiE vs. bellenguez w. all annotations ----

plot_corr <- function (df1, df2, col, title) {
  for (i in c(1,3,5,7)){
    plot(df1[df1$SNP==i,]$PRS, df2[df2$SNP == i,]$PRS, cex = 0.1, xlab = "functional finemapping", ylab = "SuSiE", 
         main = c("max SNP:",i), cex.main=1)
    abline(a=0, b=1,col = "red")
  }
        
  mtext(title,                   # Add main title
        side = 3,
        line = -1.5,
        outer = TRUE)

}


par(mfrow=c(2,2))

plot_corr(polypred,susie_polypred, SNP, "PRS")




sort_polypred1 = subset(percentage_subset(polypred1,10), percentage_subset(polypred1,10)$Diagnosis != -1)
sort_polypred3 = subset(percentage_subset(polypred3,10), percentage_subset(polypred3,10)$Diagnosis != -1)
sort_polypred5 = subset(percentage_subset(polypred5,10), percentage_subset(polypred5,10)$Diagnosis != -1)
sort_polypred7 = subset(percentage_subset(polypred7,10), percentage_subset(polypred7,10)$Diagnosis != -1)

sort_susie_polypred1 = subset(percentage_subset(susie_polypred1,10), percentage_subset(susie_polypred1,10)$Diagnosis != -1)
sort_susie_polypred3 = subset(percentage_subset(susie_polypred3,10), percentage_subset(susie_polypred3,10)$Diagnosis != -1)
sort_susie_polypred5 = subset(percentage_subset(susie_polypred5,10), percentage_subset(susie_polypred5,10)$Diagnosis != -1)
sort_susie_polypred7 = subset(percentage_subset(susie_polypred7,10), percentage_subset(susie_polypred7,10)$Diagnosis != -1)


par(mfrow=c(2,2))
plot(sort_polypred1$PRS, sort_susie_polypred1$PRS,cex = 0.1, xlab = "functional finemapping", ylab = "SuSiE", main  = "max SNP = 1")
plot(sort_polypred3$PRS, sort_susie_polypred3$PRS,cex = 0.1, xlab = "functional finemapping", ylab = "SuSiE", main  = "max SNP = 3")
plot(sort_polypred5$PRS, sort_susie_polypred5$PRS,cex = 0.1, xlab = "functional finemapping", ylab = "SuSiE", main  = "max SNP = 5")
plot(sort_polypred7$PRS, sort_susie_polypred7$PRS,cex = 0.1, xlab = "functional finemapping", ylab = "SuSiE", main  = "max SNP = 7")

# mtext("PRS ( > 2 SD)",                   # Add main title
#       side = 3,
#       line = -1.5,
#       outer = TRUE)

plot_corr(polypred_no_NA_Diagnosis, susie_polypred_no_NA_Diagnosis, PRS)



par(mfrow=c(1,1))
plot(sort_polypred1$PRS, sort_susie_polypred1$PRS,col = factor(sort_susie_polypred1$Diagnosis), cex = 0.3, xlab = "functional finemapping", ylab = "SuSiE", main  = "max SNP = 1")
 
abline(a=0, b=1,col = "red")



similarity <-lm(polypred[polypred$SNP == 3,]$PRS ~ polypred_susie[polypred_susie$SNP == 3,]$PRS)
summary(similarity)

summary(polypred[polypred$SNP == 3,]$PRS)
summary(polypred_susie[polypred_susie$SNP == 3,]$PRS)


par(mfrow=c(2,1))
plot_correlation <- function (df,main){
  plot(df[df$SNP == 1,]$PRS,  df[df$SNP == 10,]$PRS, col= ifelse(df$Diagnosis == 1, "red","black"), cex = 0.15,xlab = "max SNP = 1", ylab = "max SNP = 10", main = main)
  legend("right", legend=c("Case", "Control"),col=c("red", "black"), cex=0.8,text.font=4,inset=c(0.04,-0.0003),fill=c("red",'black'))
}
  

plot_correlation(polypred, "PolyFun")
plot_correlation(polypred_susie, "SuSiE")

## test with subset

par(mfrow=c(2,1))
plot(polypred_sd$PRS,  polypred_sd$Diagnosis, col= ifelse(polypred_sd$Diagnosis == 1, "red","black"), cex = 0.3,xlab = "PRS(PolyFun > 2SD)",ylab = "Diagnosis")
plot(polypred_susie_sd$PRS,  polypred_susie_sd$Diagnosis, col= ifelse(polypred_sd$Diagnosis == 1, "red","black"), cex = 0.3,xlab = "PRS(SuSiE >2SD)", ylab = "Diagnosis")


par(mfrow=c(2,1))
plot(polypred_sd[polypred_sd$SNP == 1,]$PRS,  polypred_susie_sd[polypred$SNP == 10,]$PRS, col= ifelse(Diagnosis == 1, "red","black"), cex = 0.15,xlab = "max SNP = 1", ylab = "max SNP = 10", main = "PolyFun")
legend("right", legend=c("Case", "Control"),col=c("red", "black"), cex=0.8,text.font=4,inset=c(0.04,-0.0003),fill=c("red",'black'))
plot(polypred_susie_sd[polypred_susie_sd$SNP == 1,]$PRS,  polypred_susie_sd[polypred_susie_sd$SNP == 10,]$PRS, col= ifelse(Diagnosis == 1, "red","black"), cex = 0.15,xlab = "max SNP = 1", ylab = "max SNP = 10", main ="SuSiE")
legend("right", legend=c("Case", "Control"),col=c("red", "black"), cex=0.8,text.font=4,inset=c(0.04,-0.0003),fill=c("red",'black'))




legend("topleft", legend=c("Line 1", "Line 2"),
       col=c("red", "blue"), lty=1:2, cex=0.8, inset=.05,fill=c("red",'black'))

# boxplot ------
par(oma = c(0,0,0,0))
par(mfrow=c(2,1))
boxplot(polypred$PRS ~ polypred$SNP, horizontal=TRUE, xlab="PRS", ylab = "num of SNP per locus", col="dodgerblue", main = "Bellenguez with PolyFun") 
#boxplot(susie_polypred$PRS ~ susie_polypred$SNP, horizontal=TRUE, xlab="PRS", ylab = "num of SNP per locus", col="dodgerblue", main = "Bellenguez with SuSiE") 
boxplot(beta_polypred$PRS ~ beta_polypred$SNP, horizontal=TRUE, xlab="PRS", ylab = "num of SNP per locus", col="dodgerblue", main = "Adjusted Beta, Polyfun") 



plot(polypred$PRS, polypred$Diagnosis,ylim = c(0,1), xlab="PRS",ylab="Diagnosis")




## multi linear regression -------


# previous try ----
# remove SNP = 1
polypred = subset(polypred[polypred$SNP !=1,])
polypred_lm_confounder <- lm(Diagnosis ~ Sex + Age, data= polypred)
summary(polypred_lm_confounder)

polypred_lm_all <- lm(Diagnosis ~ PRS + Sex + Age, data= polypred)
summary(polypred_lm_all)

polypred_lm_sd <- lm(Diagnosis ~ PRS + Sex + Age, data= polypred_sd)
summary(polypred_lm_sd) ## PRS not significant

##susie
polypred_susie_confunder <- lm(Diagnosis ~ Sex + Age , data= polypred_susie)
summary(polypred_susie_confunder)

polypred_susie_lm_all <- lm(Diagnosis ~ PRS + Sex + Age , data= polypred_susie)
summary(polypred_susie_lm_all)

### you should see each SNP seperately:----

regress_anova <- function(df){
  df = subset(df, df$Diagnosis != -1)
  polypred_lm_confounder <- lm(Diagnosis ~ Sex + Age, data= df)
  print(summary(polypred_lm_confounder))
  
  polypred_lm_all <- lm(Diagnosis ~ PRS + Sex + Age, data= df)
  print(summary(polypred_lm_all))
  
  print(anova(polypred_lm_confounder, polypred_lm_all))
}
  
regress_anova(polypred1)
regress_anova(polypred10)

regress_anova(beta_polypred1)
regress_anova(beta_polypred10)


plot(density(polypred[polypred$Diagnosis==1,]$Age),ylim=(c(0.00,0.06)), main= "Age density plot")
lines(density(polypred[polypred$Diagnosis==0,]$Age),col = "red")

polypred_susie_lm_ld<- lm(Diagnosis ~ PRS + Sex + Age, data= polypred_susie_sd)
summary(polypred_susie_lm_ld) ## PRS not significant

par(mfrow=c(1,1))
polypred10_rm = subset(polypred10,polypred10$Diagnosis!= -1)
plot(density(polypred10[polypred10$Diagnosis == 1,]$Age), col = "red", main = 'Age distribution', xlab= '', ylim = c(0,0.05))
lines(density(polypred10[polypred10$Diagnosis == 0,]$Age), main = "Control")
legend('topleft', legend=c("Case", "Control"), fill=c("Red",'Black'))


polypred10_rm$Sex = as.factor(polypred10_rm$Sex)
polypred10_rm$Diagnosis = as.factor(polypred10_rm$Diagnosis)
ggplot(polypred10_rm, aes(fill=Sex,  x=Diagnosis)) + 
  geom_bar(position="fill", stat="count")+
  ylab("percentage")+scale_y_continuous(labels = scales::percent, limits=c(0,1))+
  scale_fill_manual(values = c("lightblue","pink"),labels=c("Male","Female"))+theme_bw()



# new with filtered PRS ----

## check first
colSums(updatePRS10 == -100) ## NAs only exist in   Race and Ethnicity
summary(updatePRS10$Ethnicity)

#plot(updatePRS10$Ethnicity, xlab = "count", ylab = "Ethnicity", col=ethnicityGroup)


## running regression ----

log_reg <- function(df, population, population_subset = TRUE){
  if(population_subset == TRUE){
    df = df[df$final_population == population,]
  }
  
  mod1 <- glm(Diagnosis ~ Sex + Age, data= df, family=binomial)
  mod2 <- glm(Diagnosis ~ Sex + Age + PRS, data= df, family=binomial)
  mod3 <- glm(Diagnosis ~ Sex+ Age + X1+X2+X3+X4+X5+X6+X7+X8+X9+X10+X11+X12+X13+X14+X15+X16+X17+X18+X19+X20, data  = df, family= binomial)
  mod4 <- glm(Diagnosis ~ Sex+ Age +PRS+ X1+X2+X3+X4+X5+X6+X7+X8+X9+X10+X11+X12+X13+X14+X15+X16+X17+X18+X19+X20, data = df, family= binomial)
  
  log_result <- df %>%
    select("Sex","Age","PRS","Diagnosis")
  
  log_result <- log_result %>%
    add_column('mod1_res' = residuals(mod1, type = "pearson")) %>% 
    add_column('mod2_res' = residuals(mod2, type = "pearson")) %>%
    add_column('mod3_res' = residuals(mod3, type = "pearson")) %>% 
    add_column('mod4_res' = residuals(mod4, type = "pearson")) 
  log_result$Diagnosis<- as.factor(log_result$Diagnosis)    # Changed this into factor so that it can be matched. 
  
  par(oma = c(0,0,0,0))
  par(mfrow=c(2,1))
  
  plot(log_result$mod2_res ~ log_result$PRS, ylab = "Residual", xlab = "PRS", col = log_result$Diagnosis,
       cex=0.5, main= paste("Diagnosis ~ Sex + Age + PRS, in", population, "population"))
  
  plot(log_result$mod4_res ~ log_result$PRS, ylab = "Residual", xlab = "PRS", col = log_result$Diagnosis,
       cex=0.5, main= paste("Diagnosis ~ Sex + Age + PRS + PCs (1 to 20), in ", population, "population"))
  
  
  plot(log_result$mod2_res ~ log_result$mod1_res, ylab = "w. PRS", xlab = "w/o PRS", col = log_result$Diagnosis,
       cex=0.5,main = paste("Diagnosis ~ Sex + Age (+ PRS), in", population, "population"))
  abline(0,1, col = "blue")
  plot(log_result$mod4_res ~ log_result$mod3_res, ylab = "w. PRS", xlab = "w/o PRS", col = log_result$Diagnosis,
       cex = 0.5,main = paste("Diagnosis ~ Sex + Age (+ PRS) + PCs (1 to 20), in ", population, "population"))
  abline(0,1, col = "blue")
  
  print(summary(mod1))
  print(summary(mod2))
  print(summary(mod3))
  print(summary(mod4))
}

log_reg(updatePRS10, "EUR")
log_reg(updatePRS10, "all", population_subset = FALSE)


## strip plot -----
stripchart(SNP~PRS,
           data=polypred,
           main="PRS of different SNP per locus",
           xlab="PRS",
           ylab="SNP per locus",
           vertical=TRUE,
           pch=16
)


## violin plot -----

polypred$SNP <- as.factor(polypred$SNP)
polypred$Diagnosis_type = "Case"

polypred$Diagnosis_type[which(polypred$Diagnosis == 0)] = "Control"
polypred$Diagnosis <- as.factor(polypred$Diagnosis)

polypred_no_NA_Diagnosis <- subset(polypred, polypred$Diagnosis!= -1) 

p <- ggplot(polypred, aes(x=SNP, y=PRS)) + 
  geom_violin()


credible10 <- 
aggregrate10[c("SNP","POS",'CREDIBLE_SET')]


p + geom_boxplot(width=0.1, color="#8B4C26")+theme_bw() + labs(title = "Bellenguez + all_annotations")+
  theme(
    plot.title = element_text(hjust = 0.5, size = 14))+ xlab("Max SNP per locus")



p <- ggplot(polypred, aes(x=SNP, y=PRS, fill=Diagnosis_type))+
  theme_bw()+
  geom_violin() + labs(title = "Bellenguez + all_annotations")+
  theme(
    plot.title = element_text(hjust = 0.5, size = 14))+ xlab("Max SNP per locus")

 
par(mfrow=c(1,1))

polypred1_sd = polypred1[abs(scale(polypred1$PRS))>2,]
polypred3_sd = polypred3[abs(scale(polypred3$PRS))>2,]
polypred5_sd = polypred5[abs(scale(polypred5$PRS))>2,]
polypred7_sd = polypred7[abs(scale(polypred7$PRS))>2,]
polypred10_sd = polypred10[abs(scale(polypred10$PRS))>2,]

polypred1_sd$SNP = 1
polypred3_sd$SNP = 3
polypred5_sd$SNP = 5
polypred7_sd$SNP = 7
polypred10_sd$SNP = 10

polypred_sd = rbind(polypred1_sd,polypred3_sd,polypred5_sd,polypred7_sd,polypred10_sd)
print(dim(polypred_sd))
polypred_sd = subset(polypred_sd, polypred_sd$Diagnosis!= -1) 


polypred_sd$SNP <- as.factor(polypred_sd$SNP)
polypred_sd$Diagnosis_type = "Case"

polypred_sd$Diagnosis_type[which(polypred_sd$Diagnosis == 0)] = "Control"

p <- ggplot(polypred_sd, aes(x=SNP, y=PRS, fill=Diagnosis_type))+
  theme_bw()+
  geom_violin() + labs(title = "Bellenguez + all_annotations (> 2 SD)")+
  theme(
    plot.title = element_text(hjust = 0.5, size = 14))+ xlab("Max SNP per locus")



p<-ggplot(polypred_sd, aes(x=SNP, y=PRS, color=Diagnosis_type)) +
  geom_jitter(position=position_jitter(0.35))+theme_bw()
  labs(title = "Bellenguez + all_annotations (> 2 SD)")+
  theme(
    plot.title = element_text(hjust = 0.5, size = 14))+ xlab("Max SNP per locus")

p


stripchart(PRS~SNP,data =polypred_sd,
           main="PRS > 2SD",
           xlab="Max SNP per locus",group.names = c(1,3,5,7,10),vertical= TRUE,
           ylab="PRS",pch=16)





subset(PIP_0.5,PIP_0.5$POS %in%  PIP_0.95$POS)


aggregrate10[aggregrate10$PIP>0.95,c("SNP","BP","PIP","POS","CREDIBLE_SET")]


# comparison between plink update----
par(oma = c(0,0,4,0))
par(mfrow=c(1,3))

#max snp = 10
boxplot(PRS~Diagnosis, data = polypred10[polypred10$Diagnosis >-1,],main="previous", names=c("Control", "Case"),col = c("skyblue","lightpink"))
boxplot(PRS~Diagnosis, data = updatePRS10[updatePRS10$Diagnosis >-1,],main="updated",names =c("Control", "Case"), col= c("skyblue","lightpink"))
mtext("max snp per locus = 10", outer = TRUE, cex = 2, side = 1, line=-35)

#boxplot(PRS, data = polypred10[polypred$Diagnosis == 1,],main="previous")
#boxplot(updatePRS10$PRS,main ="updated")


#max snp= 1
boxplot(PRS~Diagnosis, data = polypred1[polypred1$Diagnosis >-1,],main="previous", names=c("Control", "Case"),col = c("skyblue","lightpink"))
boxplot(PRS~Diagnosis, data = updatePRS1[updatePRS1$Diagnosis >-1,],main="updated",names =c("Control", "Case"), col= c("skyblue","lightpink"))

mtext("max snp per locus = 1", outer = TRUE, cex = 2, side = 1, line=-35)

updatePRS10 = read.table("/gpfs/commons/home/tlin/output/prs/bellenguez_all_2/with_PC/UPDATEprs_diagnosis_0219.2021_max_snp_10_subset.tsv",header=T,fill = T)

boxplot(PRS~Diagnosis, data = polypred1[polypred1$Diagnosis >-1,],main="previous", names=c("Control", "Case"),col = c("skyblue","lightpink"),ylim = c(-0.003,0.002),cex.main=1.5)
boxplot(PRS~Diagnosis, data = updatePRS10[updatePRS10$Diagnosis!=-1,],main="PLINK updated", names=c("Control", "Case"),col = c("skyblue","lightpink"),ylim = c(-0.002,0.002), cex.main=1.5)
boxplot(PRS~Diagnosis, data = updatePRS10filter,main="remove age < 65",names =c("Control", "Case"), col= c("skyblue","lightpink"), ylim = c(-0.002,0.002),cex.main=1.5)

mtext("max snp per locus = 10", outer = TRUE, cex = 1.5, side = 1, line=-55)



## compare with sbayesR
sbayesR_bell <- read.table("/gpfs/commons/home/tlin/output/prs/sbayesr_polyfun/prs_pheno.tsv",header=T,fill = T)
sbayesR_bell

