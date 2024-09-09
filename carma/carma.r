library(CARMA)
library(data.table)
library(magrittr)
library(dplyr)
library(devtools)
library(R.utils)
library(Matrix)

args <- commandArgs(trailingOnly = TRUE)

chrom = args[1]
ld = args[2]

sumstat_path = paste("/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/carma/chr",chrom,'.tsv.gz', sep = '')
ld_path = paste('/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/LD/LD_CARMA/chr', chrom, '_', ld,'.ld',sep = '' )
output_name = paste('/gpfs/commons/home/tlin/output/CARMA/',  chrom, '_', ld, '.txt.gz', sep = '')
snp = fread(file = paste('/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/LD/LD_CARMA/chr', chrom, '_', ld,'.bim',sep = '' ), sep = "\t", select = 2)[[1]] 
print(sprintf('run CARMA using on chr %s, ld blk %s', chrom, ld))

## read data
ld = fread(file = ld_path, sep = " ", header = F, check.names = F, data.table = F, stringsAsFactors = F)
print(sprintf('%d snps in LD',dim(ld)[1]))

sumstat <- fread(file = sumstat_path , sep = "\t", header = T, check.names = F, data.table = F, stringsAsFactors = F)
sumstat <- subset(sumstat, SNP %in% snp) ## can only take the sumstat snps that's in ld blk
print(sprintf('%d snps in sumstat that overlap with LD',dim(sumstat)[1] ))

z.list<-list()
ld.list<-list()
lambda.list<-list()
z.list[[1]]<-sumstat$Z
ld.list[[1]]<-as.matrix(ld)
lambda.list[[1]]<-1

## run carma
CARMA.results<-CARMA(z.list,ld.list,lambda.list=lambda.list,outlier.switch=T)
###### Posterior inclusion probability (PIP) and credible set (CS)
sumstat.result = sumstat %>% mutate(PIP = CARMA.results[[1]]$PIPs, CS = 0)
if(length(CARMA.results[[1]]$`Credible set`[[2]])!=0){
  for(l in 1:length(CARMA.results[[1]]$`Credible set`[[2]])){ sumstat.result$CS[CARMA.results[[1]]$`Credible set`[[2]][[l]]]=l
  } }
###### write the GWAS summary statistics with PIP and CS
fwrite(x = sumstat.result,
       file = output_name, sep = "\t", quote = F, na = "NA", row.names = F, col.names = T, compress = "gzip")
