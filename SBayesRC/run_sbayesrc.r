library(SBayesRC)




ma_file="/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/sbayesrc/bellenguez_hg19.cojo"     # GWAS summary in COJO format (the only input)
ld_folder="/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/LD_sbayesrc/ukbEUR_Imputed/"        # LD reference (download from "Resources")
#annot="/gpfs/commons/home/tlin/data/sbayesrc/annot_baseline2.2.txt"         # Functional annotation (download from "Resources")
annot="/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/4SBayesRC/chr" 
suffix="_fix.tsv"
annot="$annot$chr$suffix"
echo $annot
out_prefix="/gpfs/commons/home/tlin/output/sbayesRC/bellenguez/binary_anno$chr"   # Output prefix, e.g. "./test"

##############################################
## "log2file=TRUE" means the messages will be redirected to a log file 
# Rscript -e "SBayesRC::tidy(mafile='$ma_file', LDdir='$ld_folder', \
#                   output='${out_prefix}_tidy.ma', log2file=TRUE)"
## Best practice: read the log to check issues in your GWAS summary data.   

# # Impute: optional step if your summary data doesn't cover the SNP panel
# Rscript -e "SBayesRC::impute(mafile='${out_prefix}_tidy.ma', LDdir='$ld_folder', \
#                   output='${out_prefix}_imp.ma', log2file=TRUE)"

# # SBayesRC: main function for SBayesRC
# Rscript -e "SBayesRC::sbayesrc(mafile='${out_prefix}_imp.ma', LDdir='$ld_folder', \
#                   outPrefix='${out_prefix}_imp.ma', annot='$annot', log2file=TRUE)"

SBayesRC:sbayesrc(mafile='/gpfs/commons/home/tlin/output/sbayesRC/21_imp.ma', LDdir=ld_folder, 
                 outPrefix='/gpfs/commons/home/tlin/output/sbayesRC/bellenguez/test_chr21.sbrc', annot='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/4SBayesRC/chr21_fix.tsv',
                 log2file=FALSE, bTune=FALSE)

sbayesrc = function(mafile, LDdir, outPrefix, annot="", log2file=FALSE,
                    bTune=TRUE, tuneIter=150, tuneBurn=100, thresh=0.995, 
                    tuneStep=c(0.995, 0.99, 0.95, 0.9), bTunePrior=FALSE, 
                    niter=3000, burn=1000, starth2=0.5,
                    startPi=c(0.990, 0.005, 0.003, 0.001, 0.001), 
                    gamma=c(0, 0.001, 0.01, 0.1, 1), 
                    method="sbr_ori", sSamVe="allMixVe", twopq="nbsq",
                    bOutDetail=FALSE, resam_thresh=1.1, seed=22, 
                    outFreq=10, annoSigmaScale=1.0, bOutBeta=FALSE, exclude='')