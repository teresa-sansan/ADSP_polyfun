library(coloc)
library(data.table)
install.packages(c('R.oo','R.utils'))
library(R.utils)
library(dplyr)
if(!require("remotes"))
    install.packages("remotes") # if necessary
library(remotes)
library(readr)
install.packages("gridExtra")
library(gridExtra)

install_github("chr1swallace/coloc@main",build_vignettes=TRUE)

remove_index0 = '/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/ADSP_reference_panel/fine_mapping/annotations_dl/aggregate_finemap/remove_index0/'
susie = read.table(paste(remove_index0 ,'bellenguez_susie.txt', sep = ''), header = T, sep = '\t')
dl_chr19=read.table(paste(remove_index0 ,'bellenguez_omics_dl_chr19.txt',sep=''), header = T, sep = '\t')
eqtl_microglia = fread('/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/eQTL/microglia_GENCODE_expression_full_assoc.tsv.gz', header = T, sep = '\t')
eqtl_chr19 = eqtl_microglia[eqtl_microglia$chr== 'chr19']
my.res <- coloc.abf(dataset1=susie, dataset2=eqtl_microglia)


## try their example----
data(coloc_test_data)
attach(coloc_test_data) 
par(mfrow=c(2,1))
plot_dataset(D3, main="Dataset D3")
plot_dataset(D4, main="Dataset D4")
check_dataset(D3,req="LD")
S3=runsusie(D3)


summary(test_susie)



### susieR version ----

remove.packages("susieR")

# Install the desired version from GitHub
# You may need to install the remotes package first if it's not already installed
if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
}

# Install susieR version 0.11.92 from GitHub
remotes::install_github("stephenslab/susieR@v0.11.92")

# Verify the installation
library(susieR)
packageVersion("susieR")


install.packages('susieR')
res=annotate_susie(res, snp, LD)

## reload the function from COLOC ------
annotate_susie=function(res,snp, LD) {
  ## if(ncol(res$lbf_variable) != length(snp)+1)
  ##     stop("length of snp vector should be 1 less than ncol(res$lbf_variable)")
  colnames(res$lbf_variable) = c(snp,"null")[1:ncol(res$lbf_variable)]
  colnames(res$alpha) = c(snp,"null")[1:ncol(res$alpha)]
  names(res$pip)=snp
  if(length(res$sets$cs))
    res$sets$cs = lapply(res$sets$cs, function(x) { names(x) = snp[x]; x })
  res$sld=.susie_setld(res$sets$cs,LD)
  res$pruned=FALSE
  res
}

.susie_prune=function(res,r2.prune) {
  s=res$sets$cs
  if(length(s)>1) {
    wh=which(res$sld>r2.prune,arr.ind=TRUE)
    if(length(wh)) {
      message("pruning sets in high LD")
      res = .susie_dropsets(res, wh)
      res$pruned=TRUE
    }
  }
  res
}

.susie_setld=function(s,ld) {
  ## stmp=lapply(s, setdiff, ncol(ld)+1)
  stmp=lapply(s, names)
  if(!length(stmp))
    return(0)
  if(length(stmp)==1)
    return(matrix(0,1,1,dimnames=list(names(stmp),names(stmp))))
  sld=matrix(0,length(stmp),length(stmp),dimnames=list(names(stmp),names(stmp)))
  for(i in 2:length(stmp)) {
    for(j in 1:(i-1)) {
      if(length(stmp[[i]]) && length(stmp[[j]])) {
        sld[i,j]=max(ld[stmp[[i]],stmp[[j]]]^2, na.rm=TRUE)
      } else {
        sld[i,j]=0
      }
    }
  }
  sld
}

.susie_dropsets=function(res, wh) {
  ## FIX BUG : cs_index and names of cs don't match alpha after drop - just drop from the sets
  ## drop=res$sets$cs_index[as.vector(wh)] # convert to index
  ## for(v in c("alpha","mu","mu2"))
  ##   res[[v]]=res[[v]][-drop,,drop=FALSE]
  ## for(v in c("KL","lbf","V"))
  ##   res[[v]]=res[[v]][-drop]
  drop=as.vector(wh) # raw form
  ## if("sld" %in% names(res))
  ##   res$sld=res$sld[-drop,-drop,drop=FALSE]
  res$sets$cs=res$sets$cs[-drop]
  res$sets$cs_index=res$sets$cs_index[-drop]
  res$sets$purity=res$sets$purity[-drop,,drop=FALSE]
  return(res)
}

.funrows=function(x,m,FUN) {
  x[m[1],] = apply(x[m,,drop=FALSE],2,FUN)
  x[-m[2],,drop=FALSE]
}

.sum1=function(x) {
  pmin(1, sum(x))
}

.susie_mergesets=function(res, wh) {
  m=res$sets$cs_index[as.vector(wh)] # convert to index
  res[["alpha"]] = .funrows(res[["alpha"]], m, .sum1)
  for(v in c("mu","mu2")) {
    res[[v]] = .funrows(res[[v]], m, sum)
  }
  for(v in c("KL","lbf","V"))
    res[[v]]=res[[v]][-m[2]]
  drop=as.vector(wh) # raw form
  res$sets$cs=unlist(res$sets$cs[m])
  res$sets$cs_index=res$sets$cs_index[-m[2]]
  res$sets$purity=res$sets$purity[-drop,,drop=FALSE]
  return(res)
}


###-----


vignette("a06_SuSiE",package="coloc")

coloc_path='/gpfs/commons/home/tlin/output/36k/bellenguez/adsp_ld/susie_rss/'
susie_chr12 = readRDS(paste(coloc_path, 'susie_chr12_37371914.rds', sep = ''))
coloc_chr12 = readRDS(paste(coloc_path, 'eqtl_z2_chr12_37371914_susie.rds', sep = ''))
susie_chr12_ld = read.csv(paste(coloc_path, 'susie_chr12_37371914.ld', sep = ''),sep = '\t', header = TRUE, row.names = 1)
#susie_chr12_ld_snp <- rownames(susie_chr12_ld)
susie_chr12_ld_snp <- read.csv(paste(coloc_path,'susie_chr12_37371914.tsv', sep = ''), sep ='\t')
susie_chr12_ld_snp = susie_chr12_ld_snp$SNP



annotate = annotate_susie(susie_chr12, susie_chr12_ld_snp, susie_chr12_ld)
res=coloc.susie(susie_chr12,coloc_chr12)
