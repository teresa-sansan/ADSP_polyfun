##############################################
# Variables: need to be fixed
ma_file="/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/sbayesrc/bellenguez_hg19.cojo"     # GWAS summary in COJO format (the only input)
ld_folder="/gpfs/commons/home/tlin/data/sbayesrc/ukbEUR_HM3/"        # LD reference (download from "Resources")
annot="/gpfs/commons/home/tlin/data/sbayesrc/annot_baseline2.2.txt"         # Functional annotation (download from "Resources")
out_prefix="/gpfs/commons/home/tlin/output/sbayesRC/test"   # Output prefix, e.g. "./test"
threads=4                       # Number of CPU cores

##############################################
# Code: usually don't need a change in this section
## Note: Flags were documented in the package, use ?function in R to lookup.
## We suggest to run those in multiple jobs (tasks)
export OMP_NUM_THREADS=$threads # Revise the threads

# Tidy: optional step, tidy summary data
## "log2file=TRUE" means the messages will be redirected to a log file 
Rscript -e "SBayesRC::tidy(mafile='$ma_file', LDdir='$ld_folder', \
                  output='${out_prefix}_tidy.ma', log2file=TRUE)"
## Best practice: read the log to check issues in your GWAS summary data.   

# Impute: optional step if your summary data doesn't cover the SNP panel
Rscript -e "SBayesRC::impute(mafile='${out_prefix}_tidy.ma', LDdir='$ld_folder', \
                  output='${out_prefix}_imp.ma', log2file=TRUE)"

# SBayesRC: main function for SBayesRC
Rscript -e "SBayesRC::sbayesrc(mafile='${out_prefix}_imp.ma', LDdir='$ld_folder', \
                  outPrefix='${out_prefix}_sbrc', annot='$annot', log2file=TRUE)"
# Alternative run, SBayesRC without annotation (similar to SBayesR, not recommended)
# Rscript -e "SBayesRC::sbayesrc(mafile='${out_prefix}_imp.ma', LDdir='$ld_folder', \
#                  outPrefix='${out_prefix}_sbrc_noAnnot', log2file=TRUE)"

##############################################
# Polygenic risk score
## Just a toy demo to calculate the polygenic risk score
# genoPrefix="test_chr{CHR}" # {CHR} means multiple genotype file.
## If just one genotype, input the full prefix genoPrefix="test"
# genoCHR="1-22,X" ## means {CHR} expands to 1-22 and X,
## if just one genotype file, input genoCHR=""
# output="test"
# Rscript -e "SBayesRC::prs(weight='${out_prefix}_sbrc.txt', genoPrefix='$genoPrefix', \
#                    out='$output', genoCHR='$genoCHR')"
## test.score.txt is the polygenic risk score

#################################
## SBayesRC multi
## Run each ancestry: summary data and ancestry matched LD,
##    to obtain prs1 and prs2 from the SBayesRC::prs
# prs1="eur.score.txt"
# prs2="eas.score.txt"
# tuneid="tune.id" # two columns FID IID, without header
# pheno="trait.pheno" # three columns FID IID phenotype, without header, only the samples in tuneid are used 
# outPrefix="tuned_eur_eas"
# Rscript -e "SBayesRC::sbrcMulti(prs1='$prs1', prs2='$prs2', \
#             outPrefix='$outPrefix', tuneid='$tuneid', pheno='$pheno')"
## weighted PRS in tuned_eur_eas.score.txt
## Please don't forget to exclude the tuning sample to calculate the prediction accuracy