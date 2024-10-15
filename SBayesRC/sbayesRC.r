library(SBayesRC)
library(data.table)
packageVersion('SBayesRC')
ma_file="/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/sbayesrc/rerun/bellenguez_hg38_from_parquet_reversed_a1a2.cojo"  # GWAS summary in COJO format 
#genotype="/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/LD/LD_sbayesrc/plink_file/ADSP_EUR_chr{CHR}"  # genotype prefix as LD reference (PLINK format), with {CHR} to spefify multiple
genotype="/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/LD/LD_sbayesrc/plink_file/ADSP_EUR_chr{CHR}"
outDir="/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/LD/LD_sbayesrc/ADSP_hg38/rerun_LD_sbayesrc"             # Output folder that would be created automatically
threads=4                        # Number of CPU cores for eigen decomposition


#---usually don't need change bellow
genoCHR="1-22"                     # If more than 1 genotype file, input range (e.g. "1-22") here.
# refblock=""                      # Text file to define LD blocks, by default to use our GRCH37 coordination 
refblock='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/LD/LD_sbayesrc/ADSP_hg38/ref.pos'
tool="gctb"                      # Command line to run gctb for generating the full LD matrix

##############################################
# Code
## Step1 -----
# Output $outDir/ldm.info, $outDir/ld.sh, $outDir/snplist/*.snplist

step1(mafile=ma_file, genoPrefix=genotype, outDir=outDir, genoCHR=genoCHR, blockRef=refblock, log2file= TRUE)
step1 <- function(mafile, genoPrefix, outDir, genoCHR="", blockRef="", tool="gctb", log2file=FALSE){
  output = outDir
  message("Step1: prepare the script to generate the LD matrix")
  logger.begin(output, log2file)
  if(log2file){
    message("Step1: prepare the script to generate the LD matrix")
  }
  
  chrInfo = expandCHR(genoPrefix, genoCHR)
  bMultiCHR = chrInfo$bMultiCHR
  genoFlag = checkGenoFlag(chrInfo)
  if(genoFlag != " --bfile "){
    stop("Only PLINK BED format supporetd")
  }
  
  ma = mafile
  refPos = blockRef
  if(!file.exists(ma)){
    stop("ma summary data is not there: ", ma)
  }
  
  suma = fread(ma)
  if(!"SNP" %in% colnames(suma)){
    stop("SNP must exist in ma column's name")
  }else{
    message(nrow(suma), " SNPs in summary data")
  }
  
  if(refPos == ""){
    refPos = system.file("extdata", "ref4cM_v37.pos", package="SBayesRC")
  }
  
  if(!file.exists(refPos)){
    stop("Invalid refPos: ", refPos)
  }
  
  valid_poses = fread(refPos)
  if(!all(colnames(valid_poses) %in% c("Block", "Chrom", "StartBP", "EndBP"))){
    stop("refPos is not a valid posotion file with headers: Block  Chrom  StartBP  EndBP")
  }
  setnames(valid_poses, c("blk", "chr", "start", "end"))
  
  # if(dir.exists(output)){
  #     stop("The output folder: ", output, " already exists, please change it to another one, or remove the old one")
  # }
  
  dir.create(output)
  
  
  bims = list()
  idx = 1
  for(geno1 in chrInfo$genos){
    bims[[idx]] = fread(paste0(geno1, ".bim"))
    idx = idx + 1
  }
  
  bims = rbindlist(bims)
  
  message(nrow(bims), " SNPs in genotype file")
  
  com_snp = intersect(bims$V2, suma$SNP)
  
  idx1 = match(com_snp, bims$V2)
  bims_val = bims[idx1]
  
  idx2 = match(com_snp, suma$SNP)
  ma_val = suma[idx2]
  
  m = nrow(bims_val)
  if(m == 0){
    stop("No SNPs in common between summmary data and genotype")
  }
  
  if(all(bims_val$V2 == ma_val$SNP)){
    message(m, " SNPs in common")
  }else{
    stop("Strange genotype file and summary data")
  }
  
  bA1A1 = (bims_val$V5 == ma_val$A1) & (bims_val$V6 == ma_val$A2)
  bA1A2 = (bims_val$V6 == ma_val$A1) & (bims_val$V5 == ma_val$A2)
  message("correction")
  
  bAll = bA1A1 | bA1A2
  
  message(sum(bAll), " SNPs has consistent alleles")
  
  bims = bims_val[bAll]
  if(nrow(bims) == 0){
    stop(" No SNPs left after matching the alleles")
  }
  
  for(curBlock in valid_poses$blk){
    message(" matching ", curBlock)
    curPos = valid_poses[blk==curBlock]
    chr = curPos$chr
    start = curPos$start
    end = curPos$end
    bims[V1==chr & V4>=start & V4<end, blk:=curBlock]
  }
  
  
  bims3cM = bims
  bims3cM[, newBlk:=blk]
  blks = unique(bims3cM$newBlk)
  
  outdir = output
  
  snpdir = file.path(outdir, "snplist")
  dir.create(snpdir)
  for(curBlock in blks){
    curBim = bims3cM[newBlk == curBlock]
    cat(curBim$V2, file=file.path(snpdir, paste0(curBlock, ".snplist")), sep="\n")
  }
  
  valid_poses = bims3cM[, .N, .(V1, newBlk)]
  setnames(valid_poses, "V1", "chr")
  
  valid_poses[, out:=file.path(outdir, paste0("b", newBlk))]
  valid_poses[, cmd:=paste0(tool, genoFlag, genoPrefix, " --chr ", chr, " --extract ", file.path(snpdir, paste0(newBlk, ".snplist")), " --make-full-ldm --out ", out, " &> ", out, ".log")]
  
  if(bMultiCHR){
    valid_poses[, cmd:=stringi::stri_replace_all_fixed(cmd, "{CHR}", chr)]
  }
  
  cat(valid_poses$cmd, file=file.path(outdir, "ld.sh"), sep="\n")
  
  bims[, idx:=0:(nrow(bims)-1)]
  setnames(bims, "newBlk", "Block")
  info = bims[, list(Chrom=.SD[1]$V1, StartSnpIdx = head(.SD$idx, 1), StartSnpID=head(.SD$V2, 1),  EndSnpIdx=tail(.SD$idx, 1), EndSnpID=tail(.SD$V2, 1), NumSnps=.N), by="Block"]
  
  fwrite(info, file=file.path(outdir, "ldm.info"), sep="\t", quote=F, na="NA")
  
  message("Done.")
  logger.end()
  if(log2file){
    message("Done.")
  }
}
expandCHR <- function(genoPrefix, genoCHR){
  # check if multi CHR
  bMultiCHR = FALSE
  chrs = c()
  if(genoCHR != ""){
    genoCHRnosp = gsub(" ", "", genoCHR, fixed=TRUE)
    genoStrs = c(stringi::stri_split_fixed(genoCHRnosp, ",", simplify=TRUE))
    for(genoStr in genoStrs){
      chrStr = c(stringi::stri_split_fixed(genoStr, "-", simplify=TRUE))
      if(length(chrStr) == 1){
        chrs = c(chrs, chrStr)
      }else if(length(chrStr) == 2){
        startchr = as.numeric(chrStr[1])
        endchr = as.numeric(chrStr[2])
        chrs = c(chrs, startchr:endchr)
      }else{
        stop("genoCHR doesn't recognize: ", genoCHR)
      }
    }
    
    if(!grepl("\\{CHR\\}", genoPrefix)){
      stop("genoCHR has the start and end value, however there is no {CHR} in genoPrefix string")
    }
    bMultiCHR = TRUE
  }
  
  if(!bMultiCHR){
    chrs = c("")
  }
  genos = c()
  for(curCHR in chrs){
    refGeno = gsub("{CHR}", curCHR, genoPrefix, fixed=TRUE)
    genos = c(genos, refGeno)
  }
  
  
  return(list(bMultiCHR=bMultiCHR, CHRs=chrs, genos=genos))
}


checkGenoFlag <- function(genoInfo){
  geno1 = genoInfo$genos[1]
  bBfile = TRUE
  bexts = c(".bed", ".bim", ".fam")
  for(ext in bexts){
    if(file.exists(paste0(geno1, ext))){
      bBfile = bBfile & TRUE 
    }else{
      bBfile = bBfile & FALSE
    }
  }
  
  bPfile = TRUE
  pexts = c(".pgen", ".pvar", ".psam")
  for(ext in pexts){
    if(file.exists(paste0(geno1, ext))){
      bPfile = bPfile & TRUE 
    }else{
      bPfile = bPfile & FALSE
    }
  }
  
  genoFlag = ""
  checkExts = c()
  if(bBfile & bPfile){
    warning("Found two set of genotype, use PLINK BED format")
    genoFlag = " --bfile "
    checkExts = bexts
  }else if(bBfile){
    genoFlag = " --bfile "
    message("Found PLINK BED format")
    checkExts = bexts
  }else if(bPfile){
    genoFlag = " --pfile "
    message("Found PLINK PGEN format")
    checkExts = pexts
  }else{
    stop("None genotype file found")
  }
  
  miss = c()
  for(geno1 in genoInfo$genos){
    for(ext in checkExts){
      geno1full = paste0(geno1, ext)
      if(!file.exists(geno1full)){
        miss = c(miss, geno1full)
      }
    }
  }
  if(length(miss) != 0){
    mistr = paste0(miss, collapse="\n") 
    stop("Error: some genotype file specified doesn't exist: \n", mistr)
  }
  
  return(genoFlag)
}

## step4-------------
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
outDir="/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/LD/LD_sbayesrc/ADSP_hg38/rerun_LD_sbayesrc"

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


# SBayesRC tidy the sumary data
# License: GPL
# Author: Zhili Zheng <zhilizheng@outlook.com>


#' @title Tidy GWAS summary data
#' @usage tidy(mafile, LDdir, output)
#' @param mafile string, summary data path, COJO format
#' @param LDdir string,  path to LD folder
#' @param output string, output path
#' @param freq_thresh numeric, max difference in allele frequency
#' @param N_sd_range numeric, filter the per SNP sample size in mean +- N_sd_range * sd 
#' @param rate2pq numeric, check ratio of variance estimated by 2 methods in summary data
#' @param log2file boolean, FALSE: display message on terminal; TRUE: redirect to an output file; default value is FALSE
#' @return none, results in the specified output
#' @export
tidy_fix = function(mafile, LDdir, output, freq_thresh=0.2, N_sd_range=3, rate2pq=0.5, log2file=FALSE){
  ma_file = mafile
  ld_folder = LDdir
  message('this is the right one')
  message("Tidy the summary data")
  if(file.exists(output)){
    message("the output file ", output, " exists")
    return
  }
  
  logger.begin(output, log2file)
  if(log2file){
    message("Tidy the summary data")
  }
  
  message("Summary data: ", ma_file)
  message("LD path: ", ld_folder)
  message("Output: ", output)
  message("Freq_thresh: ", freq_thresh)
  message("N_range: ", N_sd_range)
  message("rate2pq: ", rate2pq)
  
  snpinfo = fread(file.path(ld_folder, "snp.info"))
  message(nrow(snpinfo), " SNPs in LD information")
  
  ma = fread(ma_file)
  
  ma_col_names = colnames(ma)
  val_col_names = c("SNP", "A1", "A2", "freq", "b", "se", "p", "N")
  if(!all(val_col_names %in% ma_col_names)){
    stop("The summary data is not a valid COJO format")
  }
  
  ma[, p:=as.numeric(p)]
  ma[, se:=as.numeric(se)]
  ma[, b:=as.numeric(b)]
  ma[, N:=as.numeric(N)]
  ma[, freq:=as.numeric(freq)]
  
  message(nrow(ma), " SNPs in summary data")
  
  ma_val = ma[is.finite(N) & is.finite(b) & is.finite(se) & is.finite(freq) & p>=0 & p<=1]
  message(nrow(ma_val), " valid SNPs in summary data")
  
  com_snp = intersect(snpinfo$ID, ma_val$SNP)
  message(length(com_snp), " SNPs in common with LD information")
  
  idx1 = match(com_snp, snpinfo$ID)
  snpinfo_val = snpinfo[idx1]
  
  idx2 = match(com_snp, ma_val$SNP)
  ma_val2 = ma_val[idx2]
  
  if(all(ma_val2$SNP==snpinfo_val$ID)){
    message("Consistent SNP IDs")
  }else{
    stop("Duplicated SNP name in summary data")
  }
  
  bA1A1 = (ma_val2$A1 == snpinfo_val$A1) & (ma_val2$A2 == snpinfo_val$A2)
  bA1A2 = (ma_val2$A1 == snpinfo_val$A2) & (ma_val2$A2 == snpinfo_val$A1)
  
  bAll = bA1A1 | bA1A2
  
  message(sum(bAll), " SNPs have consistent alleles (A1, A2) between the summary data and LD")
  
  ma_val2[, FRQ_ref:=snpinfo_val$A1Freq]
  ma_val2[bA1A2, FRQ_ref:=1-FRQ_ref]
  
  ma_val3 = ma_val2[bAll]
  ma_val3[, frq_diff:=FRQ_ref - freq]
  ma_val4 = ma_val3[abs(frq_diff) <= freq_thresh]
  message(nrow(ma_val4), " SNPs passed the allele frequency checking with threshold ", freq_thresh)
  
  fwrite(ma_val4[, ..val_col_names], file=paste0(output, ".full"), quote=F, sep="\t", na="NA")
  
  mean_N = mean(ma_val4$N)
  sd_N =  sd(ma_val4$N)
  
  message("Mean N = ", mean_N, ", SD = ", sd_N)
  
  up_N = mean_N + N_sd_range * sd_N
  low_N = mean_N - N_sd_range * sd_N
  
  ma_val5 = ma_val4[N>=low_N & N<=up_N]
  message(nrow(ma_val5), " SNPs have sample size within mean +- ", N_sd_range, "SD")
  message("After QC, Median sample size: ", median(ma_val5$N))
  
  ma_val5[, D:=2 * freq * (1 - freq) * N]
  ma_val5[, varps:=D*(N * se^2 + b^2)/N]
  vary = median(ma_val5$varps)
  ma_val5[, indic:=sqrt( 2 * freq * (1 - freq) / vary * (N*se^2 + b^2))]
  message(" Quantile of rate2pq:")
  print(quantile(ma_val5$indic, na.rm=TRUE))
  message("Var_y: ", vary)
  ma_val6 = ma_val5[indic > (1 - rate2pq) & indic < (1 + rate2pq)] 
  message(nrow(ma_val6), " SNPs remained after QC by rate2pq ", rate2pq)
  
  fwrite(ma_val6[, ..val_col_names], file=output, quote=F, sep="\t", na="NA")
  bWarn = FALSE
  if(nrow(ma_val6) / nrow(snpinfo) < 0.7){
    bWarn = TRUE
    warning(paste0("Too many SNPs (>30%, precisely ", nrow(ma_val6)* 100 / nrow(snpinfo),"(%, were missing in the summary data after QC. The SBayesRC results may be unreliable."))
  }
  message("Done")
  logger.end()
  if(log2file){
    if(bWarn){
      warning(paste0("Too many SNPs (>30%, precisely ", nrow(ma_val6)* 100 / nrow(snpinfo),"(%) were missing in summary data after QC. The SBayesRC results may be unreliable."))
    }
    message("Done.")
  }
}


