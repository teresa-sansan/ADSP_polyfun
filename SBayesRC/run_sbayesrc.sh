#!/bin/sh
#SBATCH --job-name=sbayesrc
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=50G
#SBATCH --time=60:00:00
#SBATCH --output=/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/sbayesrc/rerun/impute.log


###SBATCH --output=/gpfs/commons/home/tlin/output/sbayesRC/bellenguez/bl/%x_%j.log
###SBATCH --array=15-22%7

module purge   
module load R/4.3.3   
module load gcc/11.2.0

#chr=$SLURM_ARRAY_TASK_ID

#Rscript /gpfs/commons/home/tlin/script/SBayesRC/run_sbayesrc.r $chr 

Rscript /gpfs/commons/home/tlin/script/SBayesRC/run_sbayesrc.r