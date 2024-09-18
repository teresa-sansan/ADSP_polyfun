devtools::install_github("ZikunY/CARMA")
library(CARMA)

pkgs = c('data.table', 'magrittr', 'dplyr','devtools', 'R.utlis')
pkgs.na = pkgs[!pkgs %in% installed.packages()[,'Package']]
if (length(pkgs.na)>0){
  install.packages(pkgs.na)
}

library(data.table)
library(magrittr)
library(dplyr)
library(devtools)
library(R.utils)
library(Matrix)


setwd('/gpfs/commons/home/tlin/data/CARMA_try')
sumstat<- fread(file = "Sample_data/sumstats_chr1_200937832_201937832.txt.gz",
                sep = "\t", header = T, check.names = F, data.table = F,
                stringsAsFactors = F)

ld = fread(file = "Sample_data/sumstats_chr1_200937832_201937832_ld.txt.gz", sep = "\t", header = F, check.names = F, data.table = F, stringsAsFactors = F)

print(head(sumstat))

z.list<-list()
ld.list<-list()
lambda.list<-list()
z.list[[1]]<-sumstat$Z
ld.list[[1]]<-as.matrix(ld)
lambda.list[[1]]<-1 
CARMA.results<-CARMA(z.list,ld.list,lambda.list=lambda.list,
                                         outlier.switch=F)
###### Posterior inclusion probability (PIP) and credible set (CS)
sumstat.result = sumstat %>% mutate(PIP = CARMA.results[[1]]$PIPs, CS = 0) 
if(length(CARMA.results[[1]]$`Credible set`[[2]])!=0){
  for(l in 1:length(CARMA.results[[1]]$`Credible set`[[2]])){ sumstat.result$CS[CARMA.results[[1]]$`Credible set`[[2]][[l]]]=l
  } }
###### write the GWAS summary statistics with PIP and CS
fwrite(x = sumstat.result,
       file = "Sample_data/sumstats_chr1_200937832_201937832_carma.txt.gz", sep = "\t", quote = F, na = "NA", row.names = F, col.names = T, compress = "gzip")
