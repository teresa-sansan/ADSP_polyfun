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


ggplot(pT[,c(prs_var, "Diagnosis",covariant)])+
       aes(x=value) +geom_histogram() 
library(reshape2)
pT_plotvar = melt(pT[,pt_prs])
ggplot(aes(x=value, colour=variable), data=pT_plotvar)+geom_density()+
  theme_light() + theme_bw()+
  scale_color_discrete("pT threshold",labels = c("0.001","0.01","0.05","0.1","0.2")) +
  ggtitle("Clumping, Thresholding, and Plumbing")+xlab('PRS')
  




