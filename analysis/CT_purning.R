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
height <- read.delim("/gpfs/commons/home/tlin/data/plink_tutorial/height_test/only_change_PLINK_input/height_prs")

## check if the order of sample ID is identical in pT and PRS(1,3,5,7,10), 
## so that we can merge different PRS into same file for easier comparison
## (Already confirmed that they have the same number of rows (11721) after removing younger and undiagnosis data)

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
plot(density(pT$PRS_5), main = "pT0.5")


plot(density(pT_check$prs_5), main = "0.5")


plot(density(height$prs_1), main = "0.1")
plot(density(height$prs_3), main = "0.3")
plot(density(height$prs_001), main = "0.001")
plot(density(height$prs_005), main = "0.005")

plot(density(height$prs_1), main = "0.1")
plot(density(height$prs_01), main = "0.01")
plot(density(height$prs_001), main = "0.001")
plot(density(height$prs_5), main = "0.5")


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
plump_thres <- c(0.001, 0.01,0.05, 0.2,0.5)
pvalue <- c(0.63852, 0.11226, 0.10054,0.10092, 0.0091)
plot(plump_thres,pvalue, xlab = "threshold", ylab = "p-value", xlim = c(0,0.55))
text(plump_thres, pvalue,c(pvalue) , pos=4, col = "blue", cex = 0.8)

