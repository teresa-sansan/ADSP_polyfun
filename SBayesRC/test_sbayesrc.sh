#!/bin/bash
#SBATCH --job-name=sbayesc_score
#SBATCH --partition=pe2
#SBATCH --nodes=1           # minimum number of nodes to be allocated
#SBATCH --ntasks=1          # number of tasks
#SBATCH --cpus-per-task=8   # number of cores on the CPU for the task
#SBATCH --mem=50G
#SBATCH --time=19:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --output=/gpfs/commons/home/tlin/output/sbayesRC/bellenguez_whole_genome/%x_%j.log

module load R/4.3.3

####SBATCH --array=10-22%6
chr=$SLURM_ARRAY_TASK_ID
###SBATCH --array=21

#bash create_cojo.sh $chr
##############################################


# Variables: need to be fixed
#ma_file='/gpfs/commons/home/tlin/data/sbayesrc/example.ma'
ma_file="/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/sbayesrc/bellenguez_hg19_chr"     # GWAS summary in COJO format (the only input)
#ma_file="$ma_file$chr.cojo"
ma_file_whole_genome='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/sbayesrc/bellenguez_hg19.cojo'
ld_folder="/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/LD_sbayesrc/ukbEUR_Imputed/"        # LD reference (download from "Resources")
#annot="/gpfs/commons/home/tlin/data/sbayesrc/annot_baseline2.2.txt"         # Functional annotation (download from "Resources")
# annot="/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/4SBayesRC/chr" 
# suffix=".tsv"
# annot="$annot$chr$suffix"
annot_whole_genome='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/4SBayesRC/fix_all.tsv'
echo $annot
out_prefix="/gpfs/commons/home/tlin/output/sbayesRC/bellenguez_whole_genome/bellenguez"
#out_prefix="/gpfs/commons/home/tlin/output/sbayesRC/bellenguez/chr$chr"   # Output prefix, e.g. "./test"
threads=4                       # Number of CPU cores


##############################################
# Code: usually don't need a change in this section
## Note: Flags were documented in the package, use ?function in R to lookup.
## We suggest to run those in multiple jobs (tasks)
export OMP_NUM_THREADS=$threads # Revise the threads

# Tidy: optional step, tidy summary data
## "log2file=TRUE" means the messages will be redirected to a log file 
# Rscript -e "SBayesRC::tidy(mafile='$ma_file_whole_genome', LDdir='$ld_folder', \
#                   output='${out_prefix}_tidy.ma', log2file=TRUE)"
# Best practice: read the log to check issues in your GWAS summary data.
# 
#  Impute: optional step if your summary data doesn't cover the SNP panel
# Rscript -e "SBayesRC::impute(mafile='${out_prefix}_tidy.ma', LDdir='$ld_folder', \
#                   output='${out_prefix}_imp.ma', log2file=TRUE)"

# # SBayesRC: main function for SBayesRC

# Rscript -e "SBayesRC::sbayesrc(mafile='${out_prefix}_imp.ma', LDdir='$ld_folder', \
#                   outPrefix='${out_prefix}_tune.sbrc', annot='$annot_whole_genome', log2file=TRUE)"

# Rscript -e "SBayesRC::sbayesrc(mafile='${out_prefix}_imp.ma', LDdir='$ld_folder', \
#                   outPrefix='${out_prefix}_imp.ma', annot='$annot', log2file=TRUE)"

# Rscript -e "SBayesRC::sbayesrc(mafile='${out_prefix}_imp.ma', LDdir='$ld_folder', \
#                  outPrefix='${out_prefix}_test.sbrc', annot='$annot', log2file=TRUE, ,bTune=FALSE)"


# Alternative run, SBayesRC without annotation (similar to SBayesR, not recommended)
# Rscript -e "SBayesRC::sbayesrc(mafile='${out_prefix}_imp.ma', LDdir='$ld_folder', \
#                   outPrefix='${out_prefix}_sbrc_noAnnot', log2file=TRUE)"

##############################################
# Polygenic risk score
## Just a toy demo to calculate the polygenic risk score
genoPrefix='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_QC/annotated_hg37_plink_ibd/qc/qc_chr{CHR}' # {CHR} means multiple genotype file.
## If just one genotype, input the full prefix genoPrefix="test"
# genoCHR="1-22,X" ## means {CHR} expands to 1-22 and X,
genoCHR="1-22"
## if just one genotype file, input genoCHR=""
# output="test"
Rscript -e "SBayesRC::prs(weight='${out_prefix}_tune.sbrc.txt', genoPrefix='$genoPrefix', \
                   out='${out_prefix}', genoCHR='$genoCHR')"
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