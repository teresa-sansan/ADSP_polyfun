
# Run in R to install dependencies
##install.packages(c("Rcpp", "data.table", "stringi", "BH",  "RcppEigen"))
# install.packages(c("Rcpp", "data.table", "stringi", "BH",  "RcppEigen"))
# # Install SBayesRC package
# install.packages("https://github.com/zhilizheng/SBayesRC/releases/download/v0.2.6/SBayesRC_0.2.6.tar.gz",
#                  repos=NULL, type="source")
# # If R report problem when installing, try alternative version (worse performance and an old version)
# install.packages("https://github.com/zhilizheng/SBayesRC/releases/download/v0.2.0/SBayesRC_0.2.0_comp.tar.gz", repos=NULL, type="source")

# library(SBayesRC)
# packageVersion(SBayesRC)

library(SBayesRC)
library(data.table)
# for (package_name in sort(loadedNamespaces())) {
#     print(paste(package_name, packageVersion(package_name)))
# }
logger.begin <- function(outPrefix, log2file){
  bLoggerTOKEN <<- log2file
  if(bLoggerTOKEN){
    zzLOGGER <<- file(paste0(outPrefix, ".log"), "wt")
    message("The messages are redirected to ", outPrefix, ".log")
    message("  No output here...")
    sink(zzLOGGER)
    sink(zzLOGGER, type="message")
  }
}

logger.end <- function(){
  if(bLoggerTOKEN){
    sink()
    sink(type="message")
    close(zzLOGGER) 
    message("  Output reactivated.")
  }
}
outDir="/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/LD/LD_sbayesrc/ADSP_hg38/LD_sbayesrc"

LDstep4 <- function(outDir, log2file=FALSE){
    gDir = outDir
    message("Step4: merge SNP information")
    message("If the LD hasn't been generated, please wait.")
    ldm = fread(file.path(gDir, "ldm.info"), head=TRUE)
    ldm[, file:=file.path(gDir, paste0("b", Block, ".ldm.full.info"))]

    infofile = file.path(gDir, "snp.info")
    if(file.exists(infofile)){
        message("snp.info generated before, the LD seems to be completed, nothing need to be done")
        return
    }

    logger.begin(paste0(infofile), log2file)
    if(log2file){
        message("Step4: merge SNP information")
    }

    infos = list()

    for(idx in 1:nrow(ldm)){
        curInfoFile = ldm$file[idx]
        if(file.exists(curInfoFile)){
            infos[[idx]] = fread(curInfoFile)
            infos[[idx]][, Block:=idx]
        }else{
            stop("Info file: ", curInfoFile, " can't be read, please check step2")
        }
    }

    info = rbindlist(infos)
    info[, Index:=0:(nrow(info)-1)]

    #info_nodup = info[!duplicated(ID)]
    
    setnames(info, c("A2", "A1", "A2Freq"), c("A1", "A2", "A1Freq"))
    
    fwrite(info[, .(Chrom, ID, Index, GenPos, PhysPos, A1, A2, A1Freq, N, Block)],  file = infofile, sep = "\t")
    #fwrite(info[, .(Chrom,ID,Index,GenPos,PhysPos,A1,A2,A2Freq,N, Block)], file=infofile, sep="\t")
    message("snp.info has been generated")
    message("To save the storage, you could remove the ", file.path(gDir, "*.ldm.bin"), " after everything cheked OK")
    message("The LD is completed!")
    logger.end()
    if(log2file){
        message("The LD is completed!")
    }
}

LDstep4(outDir=outDir, log2file=TRUE)
