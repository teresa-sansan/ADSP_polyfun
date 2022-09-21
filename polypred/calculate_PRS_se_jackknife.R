## Purpose: Calculate the mean and SE for the Jackknife PRS (replace bootstraping.)
library(readr)


calculate_se_mean <- function(polypred_path){
  print(paste('calculate PRS mean and SE from', polypred_path))
  polypred_file <- read.csv(polypred_path, header=T, sep = '\t')
  polypred_file$se <- apply(polypred_file[,c(-1,-2)], 1, sd)
  polypred_file$mean <- apply(polypred_file[,c(-1,-2)], 1, mean)
  print(dim(polypred_file))
  write.table(polypred_file[, c(1,2,dim(polypred_file)[2]-1,dim(polypred_file)[2])], file=str_replace(polypred_path,'jk','Mean_SE'), quote = F, sep = '\t', row.names = F)
}


## import polypred file
#polypred_path <- commandArgs(trailingOnly = TRUE)
all_polypred_path = c('/gpfs/commons/home/tlin/output/kunkle/kunkle_fixed_0224/finemap/polypred/10.prs_jk',
                  '/gpfs/commons/home/tlin/output/bellenguez/bellenguez_fixed_0224/finemap/max_snp_10/try_rescue_not_converge/polypred/polypred_genomewide.tsv.prs_jk',
                  '/gpfs/commons/home/tlin/output/wightman/wightman_fixed_0224/finemap/polypred/max_snp_10_new_beta_polypred.tsv.prs_jk',
                  '/gpfs/commons/home/tlin/output/jansen/finemap/polypred/max_snp_10_polypred.tsv.prs_jk')


for (polypred_path in all_polypred_path){
  calculate_se_mean(polypred_path)
}

kunkle_jk = pre_process('/gpfs/commons/home/tlin/output/prs/sumstat_jk/kunkle.tsv')
bellenguez_jk = pre_process('/gpfs/commons/home/tlin/output/prs/sumstat_jk/bellenguez.tsv')
wightman_jk = pre_process('/gpfs/commons/home/tlin/output/prs/sumstat_jk/wightman.tsv')
jansen_jk = pre_process('/gpfs/commons/home/tlin/output/prs/sumstat_jk/jansen.tsv')



