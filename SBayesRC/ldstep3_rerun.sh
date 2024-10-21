#!/bin/bash
#SBATCH --job-name=step3_fixed
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=4G
#SBATCH --time=5:00:00
#SBATCH --output=/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/LD/LD_sbayesrc/ADSP_hg38/rerun_LD_sbayesrc/%x_%j.log
#SBATCH --array=7-1196%60
module load R/4.3.3
error_file=/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/LD/LD_sbayesrc/ADSP_hg38/rerun_LD_sbayesrc/block.eigen.bin_error.txt

# # Step3: eigen decomposition for each LD block
# #  Loop idx from 1 to NUM_BLOCK (591)
# #  Submit multiple jobs on your cluster / clouds instead of for loop
# #  Input depends on $outDir/ldm.info, $outDir/b$idx.ldm.full.info, $outDir/b$idx.ldm.full.bin
# #  Output $outDir/block$block.eigen.bin, $outDir/block$block.eigen.bin.log
export OMP_NUM_THREADS=$threads  # parallel computing supported in this step
idx=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $error_file)
idx=$((idx + 1))
echo idx=${idx}

ma_file="/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/sbayesrc/rerun/bellenguez_hg38_from_parquet_reversed_a1a2.cojo"     # GWAS summary in COJO format 
#genotype="/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/LD/LD_sbayesrc/plink_file/ADSP_EUR_chr{CHR}"  # genotype prefix as LD reference (PLINK format), with {CHR} to spefify multiple
genotype="/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/LD/plinkfile_hg38_rerun/ADSP_EUR_chr{CHR}"
outDir="/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/LD/LD_sbayesrc/ADSP_hg38/rerun_LD_sbayesrc"             # Output folder that would be created automatically
threads=4                        # Number of CPU cores for eigen decomposition


Rscript -e ".libPaths('/gpfs/commons/home/tlin/R/x86_64-conda-linux-gnu-library/4.3'); library(SBayesRC); SBayesRC::LDstep3(outDir='$outDir', blockIndex=$idx, log2file=TRUE)"

