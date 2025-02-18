
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
coloc.abf(dataset1=D3, dataset2=D4)
summary(test_susie)

S3=runsusie(D3)
S4=runsusie(D4)

susie.res=coloc.susie(S3,S4)


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


try_coloc.susie=function(dataset1,dataset2, back_calculate_lbf=FALSE, susie.args=list(),  ...) {
  ## if(!requireNamespace("susieR", quietly = TRUE)) {
  ##   message("please install susieR https://github.com/stephenslab/susieR")
  ##   return(NULL)
  ## }
  if("susie" %in% class(dataset1))
    s1=dataset1
  else
    s1=do.call("runsusie", c(list(d=dataset1,suffix=1),susie.args))
  if("susie" %in% class(dataset2))
    s2=dataset2
  else
    s2=do.call("runsusie", c(list(d=dataset2,suffix=2),susie.args))
  cs1=s1$sets
  cs2=s2$sets
  if(is.null(cs1$cs) || is.null(cs2$cs) || length(cs1$cs)==0 || length(cs2$cs)==0 ) {
    warning("at least one dataset has no credible sets, nothing to colocalise")
    return(data.table(nsnps=NA))
  }
  
  idx1=cs1$cs_index
  idx2=cs2$cs_index
  bf1=s1$lbf_variable[idx1,,drop=FALSE]
  bf2=s2$lbf_variable[idx2,,drop=FALSE]
  
  ret=coloc.bf_bf(bf1,bf2,...)
  ## renumber index to match
  
  # ret$summary[,idx1:=cs1$cs_index[idx1]]
  # ret$summary[,idx2:=cs2$cs_index[idx2]]
  # ret
}

logbf_to_pp=function(bf,pi, last_is_null) {
  n=if(last_is_null) {
    ncol(bf)-1 # number of snps, null is at n+1
  } else {
    n=ncol(bf) # number of snps
  }
  if(length(pi)==1) { # scalar pi
    if(pi > 1/n)
      pi=1/n
    pi=if(last_is_null) {
      c(rep(pi,n),1-n*pi)
    } else {
      rep(pi,n)
    }
  }
  if(any(pi == 0)) { # avoid division by zero
    pi[ pi==0 ] = 1e-16
    pi=pi/sum(pi)
  }
  if(last_is_null) {
    bf=bf - bf[,ncol(bf)] ## scale bf in case needed, so log BF for null = 0
  }
  priors=matrix(log(pi),nrow(bf),ncol(bf),byrow=TRUE)
  denom=matrix(apply(bf + priors,1,logsum), nrow(bf), ncol(bf))
  exp(bf + priors - denom)
}
try_coloc.bf_bf=function(bf1,bf2, p1=1e-4, p2=1e-4, p12=5e-6, overlap.min=0.5,trim_by_posterior=TRUE) {
  if(is.vector(bf1))
    bf1=matrix(bf1,nrow=1,dimnames=list(NULL,names(bf1)))
  if(is.vector(bf2))
    bf2=matrix(bf2,nrow=1,dimnames=list(NULL,names(bf2)))
  todo <- as.data.table(expand.grid(i=1:nrow(bf1),j=1:nrow(bf2)))
  todo[,pp4:=0]
  isnps=setdiff(intersect(colnames(bf1),colnames(bf2)),
                "null")
  if(!length(isnps))
    return(data.table(nsnps=NA))
  ## scale bf in case needed
  if("null" %in% colnames(bf1))
    bf1=bf1 - matrix(bf1[,"null"],nrow(bf1),ncol(bf1))
  if("null" %in% colnames(bf2))
    bf2=bf2 - matrix(bf2[,"null"],nrow(bf2),ncol(bf2))
  
  ## check whether isnps covers the signal for each trait
  pp1=logbf_to_pp(bf1,p1, last_is_null=TRUE)
  pp2=logbf_to_pp(bf2,p2, last_is_null=TRUE)
  ph0.1=if("null" %in% colnames(pp1)) { pp1[,"null"] } else { 1 - rowSums(pp1) }
  ph0.2=if("null" %in% colnames(pp2)) { pp2[,"null"] } else { 1 - rowSums(pp2) }
  #prop1=rowSums(pp1[,c(isnps),drop=FALSE]) / rowSums(pp1[,setdiff(colnames(pp1),"null"),drop=FALSE])
  prop1=rowSums(pp1[,c(isnps),drop=FALSE]) / rowSums(pp1)
  prop2=rowSums(pp2[,c(isnps),drop=FALSE]) / rowSums(pp2[,setdiff(colnames(pp2),"null"),drop=FALSE])
  if(trim_by_posterior==TRUE) {
    drop=sapply(1:nrow(todo), function(k) {
      prop1[todo$i[k]] < overlap.min | prop2[todo$j[k]] < overlap.min
    })
    if(all(drop)) {
      warning("snp overlap too small between datasets: too few snps with high posterior in one trait represented in other")
      return(
        list(summary=cbind(
          data.table(nsnps=length(isnps),
                     hit1=colnames(pp1)[apply(pp1,1,which.max)][todo$i],
                     hit2=colnames(pp2)[apply(pp2,1,which.max)][todo$j],
                     PP.H0.abf=pmin(ph0.1[todo$i],ph0.2[todo$j]),
                     PP.H1.abf=NA, PP.H2.abf=NA, PP.H3.abf=NA, PP.H4.abf=NA),
          todo[,.(idx1=i,idx2=j)]))
      )
    }
    if(any(drop))
      todo=todo[!drop,,drop=FALSE]
  }
  ## check p12
  if(length(p12)>1) { # only proceed if didn't need to trim snps
    if(length(isnps)!=ncol(bf1) || length(isnps)!=ncol(bf2))
      stop("p12 must have length 1 or equal length to p1 and p2, but different numbers of snps between datasets")
  }
  
  ## restrict bf1/2, p1, p2, for incomplete snp overlap
  if(length(isnps)!=ncol(bf1)) {
    keep=match(isnps,colnames(bf1))
    bf1=bf1[,keep,drop=FALSE]
    if(length(p1)>1)
      p1=p1[keep]
  }
  if(length(isnps)!=ncol(bf2)) {
    keep=match(isnps,colnames(bf2))
    bf2=bf2[,keep,drop=FALSE]
    if(length(p2)>2)
      p2=p2[keep]
  }
  ## sort p12 if length(p1)>1 || length(p2)>1
  if(length(p12)==1) {
    if(length(p1)>1 || length(p2)>1) {
      p1_at_dist=min(p1)
      p2_at_dist=min(p2)
      c12=p12 / p1_at_dist / p2_at_dist
      p12=p1*p2*c12
    }
  }
  
  results <- PP <- vector("list",nrow(todo))
  ## results=lapply(1:nrow(todo), function(k) {
  for(k in 1:nrow(todo)) {
    df <- data.frame(snp=isnps, bf1=bf1[todo$i[k], ], bf2=bf2[todo$j[k], ])
    df$internal.sum.lABF <- with(df, bf1 + bf2)
    my.denom.log.abf <- logsum(df$internal.sum.lABF)
    df$SNP.PP.H4 <- exp(df$internal.sum.lABF - my.denom.log.abf)
    pp.abf <- combine.abf(df$bf1, df$bf2, p1, p2, p12, quiet=TRUE)
    PP[[k]]=df$SNP.PP.H4
    if(all(is.na(df$SNP.PP.H4))) {
      df$SNP.PP.H4=0
      pp.abf[1:5]=c(1,0,0,0,0)
    }
    common.snps <- nrow(df)
    hit1=names(which.max(bf1[todo$i[k],])) #df$snp[ which.max(abs(df$bf1)) ]
    if(is.null(hit1)) {
      hit1="-"
      pp.abf[c(1,3)]=c(0,1)
    }
    hit2=names(which.max(bf2[todo$j[k],])) #df$snp[ which.max(abs(df$bf2)) ]
    if(is.null(hit2)) {
      hit2="-"
      pp.abf[c(1,2)]=c(0,1)
    }
    results[[k]]=do.call("data.frame",c(list(nsnps=common.snps, hit1=hit1, hit2=hit2),
                                        as.list(pp.abf)))
  }
  results <- as.data.table(do.call("rbind",results))
  PP <- as.data.table(do.call("cbind",PP))
  if(nrow(todo)>1)
    colnames(PP)=paste0("SNP.PP.H4.row",1:nrow(todo))
  else
    colnames(PP)="SNP.PP.H4.abf"
  
  results=cbind(results,todo[,.(idx1=i,idx2=j)])
  PP=cbind(data.table(snp=isnps),PP)
  list(summary=results,
       results=PP,
       priors=if(length(p1)==1 && length(p2)==1 && length(p12)==1) {
         c(p1=p1,p2=p2,p12=p12)
       } else {
         list(p1=p1,p2=p2,p12=p12)
       })
}
combine.abf <- function(l1, l2, p1, p2, p12, quiet=FALSE) {
  stopifnot(length(l1)==length(l2))
  lsum <- l1 + l2
  lH0.abf <- 0
  lH1.abf <- log(p1) + logsum(l1)
  lH2.abf <- log(p2) + logsum(l2)
  lH3.abf <- log(p1) + log(p2) + logdiff(logsum(l1) + logsum(l2), logsum(lsum))
  lH4.abf <- log(p12) + logsum(lsum)
  
  all.abf <- c(lH0.abf, lH1.abf, lH2.abf, lH3.abf, lH4.abf)
  my.denom.log.abf <- logsum(all.abf)
  pp.abf <- exp(all.abf - my.denom.log.abf)
  names(pp.abf) <- paste("PP.H", (1:length(pp.abf)) - 1, ".abf", sep = "")
  if(!quiet) {
    print(signif(pp.abf,3))
    print(paste("PP abf for shared variant: ", signif(pp.abf["PP.H4.abf"],3)*100 , '%', sep=''))
  }
  return(pp.abf)
}


logsum <- function(x) {
  my.max <- max(x)                              ##take out the maximum value in log form
  my.res <- my.max + log(sum(exp(x - my.max ))) 
  return(my.res)
}

logdiff <- function(x,y) {
  my.max <- max(x,y)                              ##take out the maximum value in log form
  my.res <- my.max + log(exp(x - my.max ) - exp(y-my.max))
  return(my.res)
}


### test on PICALM -----
coloc_path='/gpfs/commons/home/tlin/output/36k/bellenguez/adsp_ld/susie_rss_oct/snp_of_interest/'
baseline_PICALM = readRDS(paste(coloc_path, 'baseline_chr11_86500001_baseline.rds', sep = ''))
eqtl_PICALM = readRDS('/gpfs/commons/home/tlin/output/36k/bellenguez/adsp_ld/susie_rss_oct/eqtl_tss/eqtl_chr11_ENSG00000073921.18.rds')
eqtl_PICALM_ld = read.csv('/gpfs/commons/home/tlin/output/36k/bellenguez/adsp_ld/susie_rss_oct/eqtl_tss/eqtl_chr11_ENSG00000073921.18.ld')
PICALM_ld = read.csv(paste(coloc_path, 'baseline_chr11_84957209_baseline.ld', sep = ''),sep = '\t', header = TRUE, row.names = 1)
omics_PICALM = readRDS(paste(coloc_path, 'omics_chr11_84957209_omics.rds', sep = ''))

eqtl_PICALM_snp <- read.csv('/gpfs/commons/home/tlin/output/36k/bellenguez/adsp_ld/susie_rss_oct/eqtl_tss/eqtl_chr11_ENSG00000073921.18.tsv', sep = '\t')
eqtl_PICALM_snp = eqtl_PICALM_snp$SNP
annotate_PICALM = annotate_susie(eqtl_PICALM, eqtl_PICALM_snp, eqtl_PICALM_ld)
summary(annotate_PICALM)



bl_PICALM_snp <- read.csv(paste(coloc_path,'baseline_chr11_84957209_baseline.tsv', sep = ''), sep ='\t')
bl_PICALM_snp = bl_PICALM_snp$SNP. #11666

omics_PICALM_snp <- read.csv(paste(coloc_path,'omics_chr11_84957209_omics.tsv', sep = ''), sep ='\t')
omics_PICALM_snp = omics_PICALM_snp$SNP.  #11666

identical(omics_PICALM_snp, bl_PICALM_snp) 

annotate_bl = annotate_susie(baseline_PICALM, bl_PICALM_snp, PICALM_ld)
annotate_omics = annotate_susie(omis_PICALM, omics_PICALM_snp, PICALM_ld)

res=coloc.susie(annotate_bl,annotate_omics)
summary(annotate_bl)
summary(annotate_omics)
sensitivity(res,"H4 > 0.9",row=1,dataset1=annotate_bl,dataset2=annotate_omics)



##############




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

