#!/bin/bash
#SBATCH --job-name=step3_rerun
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=30G
#SBATCH --array=2087-2108%10
#SBATCH --time=15:00:00
#SBATCH --output=/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/LD/LD_sbayesrc/ADSP_hg38/rerun_LD_sbayesrc/%x_%j.log

##
###
# source /gpfs/commons/home/tlin/miniconda3/etc/profile.d/conda.sh
# conda activate sbayesrc
module load R/4.3.3
####SBATCH --array=11-2108%60
##############################################

#CHR=$SLURM_ARRAY_TASK_ID
ma_file="/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/sbayesrc/rerun/bellenguez_hg38_from_parquet_reversed_a1a2.cojo"     # GWAS summary in COJO format 
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
# Step1: generate the LD block information and script
# Output $outDir/ldm.info, $outDir/ld.sh, $outDir/snplist/*.snplist
# Rscript -e "SBayesRC::LDstep1(mafile='${ma_file}', genoPrefix='$genotype', outDir='$outDir', genoCHR='$genoCHR', blockRef='$refblock', log2file = TRUE)"

# Step2: generate each LD matrix for blocks
#  Loop idx from 1 to NUM_BLOCK (591)
#  Submit multiple jobs on your cluster / clouds instead of for loop
#  Input depends on $outDir/ld.sh, $outDir/snplist/$idx.snplist
#  Ouput $outDir/b$idx.ldm.full.info, $outDir/b$idx.ldm.full.bin

# idx=$SLURM_ARRAY_TASK_ID
# which gctb
#Rscript -e ".libPaths('/gpfs/commons/home/tlin/R/x86_64-conda-linux-gnu-library/4.3'); library(SBayesRC); SBayesRC::LDstep2(outDir='$outDir', blockIndex=$idx, log2file=TRUE)"


# for idx in {1..591}; do
#     # Rscript -e "SBayesRC::LDstep2(outDir='$outDir', blockIndex=$idx, log2file=TRUE)"
# done

# # Step3: eigen decomposition for each LD block
# #  Loop idx from 1 to NUM_BLOCK (591)
# #  Submit multiple jobs on your cluster / clouds instead of for loop
# #  Input depends on $outDir/ldm.info, $outDir/b$idx.ldm.full.info, $outDir/b$idx.ldm.full.bin
# #  Output $outDir/block$block.eigen.bin, $outDir/block$block.eigen.bin.log
# export OMP_NUM_THREADS=$threads  # parallel computing supported in this step
# for idx in {2087..2108}
# do
#     echo $idx
#     Rscript -e ".libPaths('/gpfs/commons/home/tlin/R/x86_64-conda-linux-gnu-library/4.3'); library(SBayesRC); SBayesRC::LDstep3(outDir='$outDir', blockIndex=$idx, log2file=TRUE)"
# done
#Rscript -e "library(SBayesRC);SBayesRC::LDstep3(outDir='$outDir', blockIndex=$idx, log2file=TRUE)"
# for idx in {1..591}; do
#     Rscript -e "SBayesRC::LDstep3(outDir='$outDir', blockIndex=$idx, log2file=TRUE)"
# done

# # Step4: merge LD information
Rscript -e ".libPaths('/gpfs/commons/home/tlin/R/x86_64-conda-linux-gnu-library/4.3'); library(SBayesRC); SBayesRC::LDstep4(outDir='$outDir', log2file=TRUE)"


# Step5: clean if necessary
# Essential for analysis: $outDir/ldm.info, $outDir/snp.info, $outDir/block*.eigen.bin 
# Other files could be removed
# Note: before removing, check all blocks were finished.