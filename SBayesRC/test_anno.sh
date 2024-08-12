#!/bin/bash
#SBATCH --job-name=sbayesr_noanno
#SBATCH --partition=pe2
#SBATCH --nodes=1           # minimum number of nodes to be allocated
#SBATCH --ntasks=1          # number of tasks
#SBATCH --cpus-per-task=8   # number of cores on the CPU for the task
#SBATCH --mem=120G
#SBATCH --time=9:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --output=/gpfs/commons/home/tlin/output/sbayesRC/bellenguez_whole_genome/%x_%j.log

module load R/4.3.3
chr=$SLURM_ARRAY_TASK_ID
###SBATCH --array=21



# Variables: need to be fixed
#ma_file='/gpfs/commons/home/tlin/data/sbayesrc/example.ma'
ma_file="/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/sbayesrc/bellenguez_hg19_chr"     # GWAS summary in COJO format (the only input)
ma_file="$ma_file$chr.cojo"
ld_folder="/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/LD_sbayesrc/ukbEUR_Imputed/"        # LD reference (download from "Resources")
#annot="/gpfs/commons/home/tlin/data/sbayesrc/annot_baseline2.2.txt"         # Functional annotation (download from "Resources")
annot="/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/4SBayesRC/new_chr" 
suffix=".tsv"
annot="$annot$chr$suffix"
echo $annot
out_prefix="/gpfs/commons/home/tlin/output/sbayesRC/bellenguez/chr$chr"   # Output prefix, e.g. "./test"
threads=4                       # Number of CPU cores


##############################################
# Code: usually don't need a change in this section
## Note: Flags were documented in the package, use ?function in R to lookup.
## We suggest to run those in multiple jobs (tasks)
export OMP_NUM_THREADS=$threads # Revise the threads

Rscript -e "SBayesRC::sbayesrc(mafile='${out_prefix}_imp.ma', LDdir='$ld_folder', \
                 outPrefix='${out_prefix}_test.sbrc', annot='$annot', log2file=TRUE, bTune=FALSE)"

