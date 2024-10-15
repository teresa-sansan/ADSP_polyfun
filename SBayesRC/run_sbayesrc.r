library(SBayesRC)
packageVersion('SBayesRC') ## make sure you are using the latest version (0.2.6)

# args <- commandArgs(trailingOnly = TRUE)
# chrom = args[1]

#ma_file=paste("/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/sbayesrc/bellenguez_hg38_chr",chrom,"_from_parquet_reversed_a1a2.cojo",sep = '')     # GWAS summary in COJO format (the only input)
ma_file="/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/sbayesrc/rerun/bellenguez_hg38_from_parquet_reversed_a1a2.cojo"
#ld_folder="/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/LD_sbayesrc/ukbEUR_Imputed/"        # LD reference (download from "Resources")
ld_folder='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/LD/LD_sbayesrc/ADSP_hg38/LD_sbayesrc'
##annot="/gpfs/commons/home/tlin/data/sbayesrc/annot_baseline2.2.txt"         # Functional annotation (download from "Resources")
#annot=paste("/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/4SBayesRC/hg38/bl_chr",chrom,".tsv", sep = '' )
#annot='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations_hg38/merged_annotations_ADSP_v2/baseline_filtered/chr22.annot.gz'
suffix=".tsv"
#annot="$annot$chr$suffix"
out_prefix="/gpfs/commons/home/tlin/output/sbayesRC/bellenguez/bl/chr$chr"   # Output prefix, e.g. "./test"

##############################################
## "log2file=TRUE" means the messages will be redirected to a log file 
## tidy sumstat
#tidy(mafile=ma_file, LDdir=ld_folder, output= paste(dirname(ma_file), '/bellenguez_hg38_chr',chrom,'_reversed_a1a2_tidy.ma', sep = ''), log2file=TRUE)
tidy(mafile=ma_file, LDdir=ld_folder, output= paste(dirname(ma_file), '/bellenguez_hg38_reversed_a1a2_tidy.ma', sep = ''), log2file=TRUE)

# # Impute: optional step if your summary data doesn't cover the SNP panel
# Rscript -e "SBayesRC::impute(mafile='${out_prefix}_tidy.ma', LDdir='$ld_folder', \
#                   output='${out_prefix}_imp.ma', log2file=TRUE)"

#impute(mafile=paste(dirname(ma_file), '/bellenguez_hg38_chr',chrom,'_reversed_a1a2_tidy.ma', sep = ''), LDdir=ld_folder,output=paste(dirname(ma_file), '/bellenguez_hg38_chr',chrom,'_reversed_a1a2_imp.ma', sep = ''), log2file=TRUE)
impute(mafile=paste(dirname(ma_file), '/bellenguez_hg38_reversed_a1a2_tidy.ma', sep = ''), LDdir=ld_folder,output=paste(dirname(ma_file), '/bellenguez_hg38_reversed_a1a2_imp.ma', sep = ''), log2file=TRUE)


# # SBayesRC: main function for SBayesRC
#sbayesrc(mafile=paste(dirname(ma_file), '/bellenguez_hg38_chr',chrom,'_reversed_a1a2_tidy.ma', sep = ''), LDdir=ld_folder, annot = annot, outPrefix=paste('/gpfs/commons/home/tlin/output/sbayesRC/bellenguez/bl/chr',chrom,'_sbrc', sep = ''), log2file=TRUE)

#sbayesrc(mafile=paste(dirname(ma_file), '/bellenguez_hg38_chr',chrom,'_reversed_a1a2_imp.ma', sep = ''), LDdir=ld_folder, annot = annot, outPrefix=paste('/gpfs/commons/home/tlin/output/sbayesRC/bellenguez/bl/chr',chrom,'_sbrc', sep = ''), log2file=TRUE)

# sbayesrc(mafile=ma_file, LDdir=ld_folder, 
#          outPrefix='/gpfs/commons/home/tlin/output/sbayesRC/hg38/test_chr22.sbrc', annot=annot,
#          log2file=FALSE, bTune=FALSE)

# SBayesRC:sbayesrc(mafile=ma_file, LDdir=ld_folder, 
#                  outPrefix='/gpfs/commons/home/tlin/output/sbayesRC/hg38/test_chr22.sbrc', annot=annot,
#                  log2file=FALSE, bTune=FALSE)

# sbayesrc = function(mafile, LDdir, outPrefix, annot="", log2file=FALSE,
#                     bTune=TRUE, tuneIter=150, tuneBurn=100, thresh=0.995, 
#                     tuneStep=c(0.995, 0.99, 0.95, 0.9), bTunePrior=FALSE, 
#                     niter=3000, burn=1000, starth2=0.5,
#                     startPi=c(0.990, 0.005, 0.003, 0.001, 0.001), 
#                     gamma=c(0, 0.001, 0.01, 0.1, 1), 
#                     method="sbr_ori", sSamVe="allMixVe", twopq="nbsq",
#                     bOutDetail=FALSE, resam_thresh=1.1, seed=22, 
#                     outFreq=10, annoSigmaScale=1.0, bOutBeta=FALSE, exclude='')
