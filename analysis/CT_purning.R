library(ggplot2)
library(lattice)
library(dplyr)
library(stringr)
library(tidyr)
library(ggstance)
library(jtools)
library(tibble)
library(readr)

pT <- pre_process("/gpfs/commons/home/tlin/output/prs/bellenguez_all_2/pT/pT_PRS_withPC.tsv")
pT_check <- read.delim('/gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap/pT/check_prs')
all(pT$SampleID == updatePRS1$SampleID)
all(pT$SampleID == updatePRS10$SampleID)

pT = pT %>%
     add_column(prs_max1 = updatePRS1$PRS,
             prs_max3 = updatePRS3$PRS,
             prs_max5 = updatePRS5$PRS,
             prs_max7 = updatePRS7$PRS,
             prs_max10 = updatePRS10$PRS)

ggplot(pT, aes(x=prs_max10)) + geom_density()


prs_var = c("prs_max1","prs_max5","prs_max10")
pt_prs= c("PRS_001","PRS_01","PRS_05","PRS_1","PRS_2")
covariant = c("Sex","APOE","Age","final_population")



par(mfrow=c(2,2))
plot(density(pT$PRS_001), main = "0.001")
plot(density(pT$PRS_01), main = "0.01")
plot(density(pT$PRS_1), main = "0.1")
plot(density(pT$PRS_5), main = "0.5")

plot(density(pT$prs_max1), main = "max1")
plot(density(pT$prs_max3), main = "max3")
plot(density(pT$PRS_1), main = "pT_0.1")
plot(density(pT$PRS_5), main = "pT_0.5")


plot(density(pT_check$prs_5), main = "0.5")
plot(density(pT_check$prs_1), main = "0.1")
plot(density(pT_check$prs_01), main = "0.01")
plot(density(pT_check$prs_001), main = "0.001")

plot(density(height$prs_1), main = "height_pt_0.1")
plot(density(height$prs_3), main = "height_pt_0.3")
plot(density(height$prs_001), main = "height_pt_0.001")
plot(density(height$prs_005), main = "height_pt_0.005")



ggplot(pT[,c(prs_var, "Diagnosis",covariant)])+
       aes(x=value) +geom_histogram() 
library(reshape2)
pT_plotvar = melt(pT[,pt_prs])
ggplot(aes(x=value, colour=variable), data=pT_plotvar)+geom_density()+
  theme_light() + theme_bw()+
  scale_color_discrete("pT threshold",labels = c("0.001","0.01","0.05","0.1","0.2")) +
  ggtitle("Clumping and Thresholding")+xlab('PRS')
  

mod <- glm(Diagnosis ~ Sex + Age + PRS_1, data= pT, family=binomial)

residuals(mod, type = "pearson")
summary(mod)

par(mfrow=c(1,1))
p_thres <- c(0.001, 0.01,0.05, 0.2,0.5)
pvalue <- c(0.63852, 0.11226, 0.10054,0.10092, 0.0091)
plot(plump_thres,pvalue, xlab = "threshold", ylab = "p-value", xlim = c(0,0.55))
text(plump_thres, pvalue,c(pvalue) , pos=4, col = "blue", cex = 0.8)

for(i in p_thres){
  print(paste("p value threshold = ",i))
  i = paste("PRS",str_replace(i,'0.','_'), sep='')


  mod <- glm(substitute(Diagnosis ~ Sex + Age +i ,list(i = as.name(i))), data= pT, family=binomial)
  #prs.beta <- as.numeric(mod$coefficients[1])
  #prs.se <- as.numeric(mod$coefficients[2])
  #prs.p <- as.numeric(mod$coefficients[4])
  #prs.result <- rbind(prs.result, data.frame(Threshold=i, P=prs.p, BETA=prs.beta,SE=prs.se))
  print(summary(mod))
}


## AUC -----

library(pROC)
roc_max1 <- roc(pT$Diagnosis, pT$prs_max1)

roc_max3 <- roc(pT$Diagnosis, pT$prs_max3)
plot(roc_max3)

roc_PT05 <-  roc(pT$Diagnosis, pT$PRS_5)
plot(roc_PT05)

rocobj1 <-plot.roc(pT$Diagnosis, pT$prs_max1, main = 'AUC comparison')
rocobj2 <-lines.roc(pT$Diagnosis, pT$prs_max3)
testobj <- roc.test(rocobj1,rocobj2)
text(50,50,labels=paste("pvalue = ", format.pval(testobj$p.value)), adj=c(0,.5))

plot.roc(roc_max1, print.auc=TRUE)
plot.roc(roc_max3, print.auc=TRUE)

 ggroc(list(max1=roc_max1, max3=roc_max3))
 
 
roc_list <- roc(Diagnosis ~prs_max1+prs_max5+prs_max10+PRS_01+PRS_5, data = pT)

cal_auc <- function(col, df = pT){
  AUC = round(auc(pT$Diagnosis, pT$col),3)
  return(AUC)
}

rocvalue<-list()

for (i in list("prs_max1","prs_max5","prs_max10","PRS_01","PRS_5")){
  print (i)
  print(substitute(auc(pT$Diagnosis, pT$i ,list(i = as.name(i)))))
  }


curve <- ggroc(roc_list)
curve + xlab("FPR") + ylab("TPR") +
  scale_colour_discrete(name="PRS",labels = c(paste("max_snp = 1, AUC = ",cal_auc(prs_max1)cal_auc(prs_max1)), "Treatment 1", "1","2","3"))
  
                        