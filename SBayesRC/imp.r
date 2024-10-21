library(SBayesRC)
packageVersion('SBayesRC') ## make sure you are using the latest version (0.2.6)

library(Rcpp)
sourceCpp("getLDPrefix.cpp")
exists("getLDPrefix")  # This should return TRUE



args <- commandArgs(trailingOnly = TRUE)
chrom = args[1]
print(paste('running on chr', chrom))

ma_file=paste("/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/sbayesrc/rerun/bellenguez_hg38_chr",chrom,"_from_parquet_reversed.cojo",sep = '')     # GWAS summary in COJO format (the only input)
#ma_file="/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/sbayesrc/rerun/bellenguez_hg38_from_parquet_reversed_a1a2.cojo"
#ld_folder="/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/LD_sbayesrc/ukbEUR_Imputed/"        # LD reference (download from "Resources")
ld_folder='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/LD/LD_sbayesrc/ADSP_hg38/rerun_LD_sbayesrc'
##annot="/gpfs/commons/home/tlin/data/sbayesrc/annot_baseline2.2.txt"         # Functional annotation (download from "Resources")
#annot=paste("/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/4SBayesRC/hg38/bl_chr",chrom,".tsv", sep = '' )
#annot='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations_hg38/merged_annotations_ADSP_v2/baseline_filtered/chr22.annot.gz'
suffix=".tsv"
#annot="$annot$chr$suffix"
#out_prefix="/gpfs/commons/home/tlin/output/sbayesRC/bellenguez/bl/chr$chr"   # Output prefix, e.g. "./test"


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

# getLDPrefix2 <- function(mldm) {
#     .Call(`_SBayesRC_getLDPrefix2`, mldm)
# }


# snpinfo = fread(file.path(ld_folder, "snp.info.all"))
# snpinfo = snpinfo[snpinfo$Chrom == chrom,]
# write.table(snpinfo, paste("/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/LD/LD_sbayesrc/ADSP_hg38/rerun_LD_sbayesrc/chr", chrom, "snp.info", sep = ''), row.names = FALSE, sep = "\t", quote = FALSE)

impute_fix = function(mafile, LDdir, output, chrom, thresh=0.995, log2file=FALSE){
    LD_folder = LDdir
    message("Impute the summary data by LD")
    if(file.exists(output)){
        message("don't need to run: ", output, " exists")
        return
    }

    logger.begin(output, log2file)
    if(log2file){
        message("Impute the summary data by LD")
    }

    ma = fread(mafile)
    ma[, D:=2*freq * (1-freq) * N]
    ma[, varps:=D*(N * se^2 + b^2)/N]

    vp_median = median(ma$varps)
    N_median = median(ma$N)


    snpinfo = fread(file.path(LD_folder, "snp.info.all"), head=TRUE)
    if(ncol(snpinfo) == 9){
        setnames(snpinfo, c("CHR", "SNP", "GD", "BP", "A1", "A2", "freq", "N", "Block"))
    }else if(ncol(snpinfo) == 8){
        setnames(snpinfo, c("CHR", "SNP", "BP", "A1", "A2", "freq", "N", "Block"))
    }else if(ncol(snpinfo) == 10){
        setnames(snpinfo, c("CHR", "SNP", "Index", "GD", "BP", "A1", "A2", "freq", "N", "Block"))
    }else{
        stop("the LD information looks odd")
    }
    snpinfo=snpinfo[snpinfo$CHR == chrom,]
    snpinfo[, .(SNP, A1, A2, freq, Block)]
    snpfinal = snpinfo[, .(SNP, A1, A2, freq, Block)]

    comSNP = intersect(ma$SNP, snpfinal$SNP)

    ma_run = ma[SNP %in% comSNP]

    idx = match(ma_run$SNP, snpfinal$SNP)
    message(length(idx), " SNPs in common between GWAS summary and LD")

    bA1A1 = (ma_run$A1 == snpfinal[idx]$A1) & (ma_run$A2 == snpfinal[idx]$A2)
    message(sum(bA1A1), " SNPs set from summary data")

    bA1A2 = (ma_run$A1 == snpfinal[idx]$A2) & (ma_run$A2 == snpfinal[idx]$A1)
    message(sum(bA1A2), " SNPs flipped alleles")

    message(sum(bA1A1 | bA1A2), " SNPs are typed SNPs")

    snpfinal[, b:=NA_real_]
    snpfinal[, se:=NA_real_]
    snpfinal[, p:=NA_real_]
    snpfinal[, N:=NA_real_]
    snpfinal[idx[bA1A1], freq:=ma_run[bA1A1]$freq]
    snpfinal[idx[bA1A1], b:=ma_run[bA1A1]$b]
    snpfinal[idx[bA1A1], se:=ma_run[bA1A1]$se]
    snpfinal[idx[bA1A1], p:=ma_run[bA1A1]$p]
    snpfinal[idx[bA1A1], N:=ma_run[bA1A1]$N]

    snpfinal[idx[bA1A2], freq:=(1.0 - ma_run[bA1A2]$freq)]
    snpfinal[idx[bA1A2], b:=(-ma_run[bA1A2]$b)]
    snpfinal[idx[bA1A2], se:=ma_run[bA1A2]$se]
    snpfinal[idx[bA1A2], p:=ma_run[bA1A2]$p]
    snpfinal[idx[bA1A2], N:=ma_run[bA1A2]$N]

    snpfinal[is.na(N), N:=N_median]


    message("Start summary imputation...")
    all_ma = list()

    idxBlocks = unique(snpfinal$Block)

    info = getLDPrefix(LD_folder)
    #info=_SBayesRC_getLDPrefix(LD_folder)

    for(idxBlk in idxBlocks){
        message("==========", idxBlk, "=========")

        ma_Block = snpfinal[Block==idxBlk]
        ma_Block[, r2:=1]

        m = nrow(ma_Block)

        idxtt = ma_Block[is.finite(b), which=TRUE]

        if(length(idxtt) != 0){
            ma_exist = ma_Block[idxtt]
            ma_exist[, z:=b/se]

            message(" Imputing... ")
            Zi = impGa(info$template, idxBlk, info$type, ma_exist$z, idxtt, m, thresh)

            ma_block_n = ma_Block[-idxtt]
            ma_block_n[, z2:=Zi]
            ma_block_n[, base1:=sqrt(2*freq*(1-freq)*(N+z2*z2))]
            ma_block_n[, b2:=z2*sqrt(vp_median)/base1]
            ma_block_n[, se2:=sqrt(vp_median)/base1]
            ma_block_n[, p2:=pchisq(z2*z2, df=1, lower.tail=FALSE)]
            ma_block_n[, r2:=0] # indicate imputed

            ma_Block[is.na(p), b:=ma_block_n$b2]
            ma_Block[is.na(p), se:=ma_block_n$se2]
            ma_Block[is.na(p), r2:=ma_block_n$r2]
            ma_Block[is.na(p), p:=ma_block_n$p2]
        }else{
            ma_Block[is.na(p), b:=0]
            ma_Block[is.na(p), se:=1]
            ma_Block[is.na(p), r2:=-1]
            ma_Block[is.na(p), p:=1]
        }


        all_ma[[idxBlk]] = ma_Block
    }

    ma_out = rbindlist(all_ma)
    ma_out[, Block:=NULL]
    fwrite(ma_out, file=output, sep="\t", quote=F, na="NA")

    message("Done.")
    print(proc.time())
    logger.end()
    if(log2file){
        message("Done.")
    }
}


impute_fix(mafile=paste(dirname(ma_file), '/bellenguez_hg38_chr',chrom,'_reversed_a1a2_tidy.ma', sep = ''), LDdir=ld_folder, chrom=chrom,
        output=paste(dirname(ma_file), '/bellenguez_hg38_chr',chrom,'_reversed_a1a2_imp.ma', sep = ''))


